import subprocess
import re
import tempfile
import csv
import glob
import json
import os
import time
import random
from collections import OrderedDict
from io import StringIO
from lxml import etree
import sys
import requests

import shutil
import urllib

def safe_read(element, xpath, index=None, xpath_fallback=None):

    # the logic here is basically, try one tag, if it's not there, return an
    # alternate tag. for either case, if there's a specific element you want
    # (like 0), just return that instead of the list.

    # also, I haven't tested to be sure this is really all that safe.  it's just
    # safer than the way I was doing it.

    val = element.xpath(xpath)
    if val:
        if index is not None:
            return val[index]
        else:
            return val
    elif xpath_fallback is not None:
        return safe_read(element, xpath_fallback, index)
    else:
        return ''


def get_accession_metadata(accession_id, sra_metadata_file):
    print("Getting accession: {}".format(str(accession_id)))
    params = { 'db': 'sra', 'rettype': 'docset', 'id': accession_id }
    # params = { 'save': 'efetch', 'db': 'sra', 'rettype': 'docset', 'term': accession_id }
    retry_count = 0
    while True:
        ret = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi', params=params)

        # ret = requests.get('https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi', params=params)
        if ret.status_code == 429:
            delay = retry_count + random.uniform(0, 2)
            print(  "Delaying for 429 error " + str(delay), file=sys.stderr)
            time.sleep(delay)
            retry_count = retry_count + 1
        elif ret.status_code != 200:
            return_code = ret.status_code
            break
        else:
            break

    #print("GOT: ", ret, params);
    parser = etree.XMLParser(remove_blank_text=True)
#    if type(ret.text) == unicode:
#        result_obj = StringIO.StringIO(ret.text.encode('utf-8'));
#    else:
#        result_obj = StringIO.StringIO(ret.text);
    tree = etree.fromstring(ret.text.encode('ascii'))

#    tree = etree.parse(result_obj, parser)

    # print it out just for fun
    if sra_metadata_file:
        fp = file(sra_metadata_file, "w")
        print( etree.tostring(tree, pretty_print=True), file=fp)
        fp.close()
    
    return parse_accession_metadata(accession_id, tree)

def parse_sra_xml(xml_fh):

    parser = etree.XMLParser(remove_blank_text=True)
    tree = etree.parse(xml_fh, parser)
    return parse_accession_metadata(None, tree)

def parse_accession_metadata(accession_id, tree):

    # step through the experiments
    exp_meta = []
    for experiment_package in tree.xpath('//EXPERIMENT_PACKAGE'):
        exp = {}
        exp['exp_id'] = safe_read(experiment_package, 'EXPERIMENT/@accession', index=0)

        plat = safe_read(experiment_package, 'EXPERIMENT/PLATFORM/*[1]', index=0)
        if plat != '':
            exp['platform_name'] = plat.tag
            exp['instrument_model'] = safe_read(plat, 'INSTRUMENT_MODEL/text()')
            if isinstance(exp['instrument_model'], list):
                exp['instrument_model'] = exp['instrument_model'][0]

        exp['design_description'] =''.join(safe_read(experiment_package, 'EXPERIMENT/DESIGN/DESIGN_DESCRIPTION/text()'))
        exp['study_id'] = safe_read(experiment_package, 'EXPERIMENT/STUDY_REF/@accession', index=0)

        for db in experiment_package.xpath('STUDY//XREF_LINK/DB/text()'):
            exp['study_'+db+'_ids'] = experiment_package.xpath('STUDY//XREF_LINK/ID/text()')

        exp['library_name'] = safe_read(experiment_package, 'EXPERIMENT//LIBRARY_NAME/text()', index=0)
        exp['library_selection'] = safe_read(experiment_package, 'EXPERIMENT//LIBRARY_SELECTION/text()', index=0)
        exp['library_strategy'] = safe_read(experiment_package, 'EXPERIMENT//LIBRARY_STRATEGY/text()', index=0)
        #this might be unreliable. use the existence of paired file
        exp['library_layout'] = "PAIRED" if safe_read(experiment_package, 'EXPERIMENT//LIBRARY_LAYOUT/PAIRED', index=0) != "" else "SINGLE" 

        exp['runs'] = []
        for run in experiment_package.xpath('RUN_SET/RUN'):
            rdata = {}
            rdata['run_id'] = safe_read(run, '@accession')[0]
            my_out = "run: {}, exp: {}, study: {}".format(rdata['run_id'], exp['exp_id'], exp['study_id'])
            if accession_id != None and rdata['run_id'] != accession_id and accession_id != exp['exp_id'] and accession_id != exp['study_id']:
                print( "Skipping -- " + my_out, file=sys.stderr)
                continue
            else:
                print( "Using -- " + my_out, file=sys.stderr)
            rdata['accession'] = rdata['run_id']
            try:
                rdata['total_bases'] = int(safe_read(run, '@total_bases')[0])
                rdata['total_spots'] = int(safe_read(run, '@total_spots')[0])
                rdata['size'] = int(safe_read(run, '@size')[0])
            except Exception as e:
                print( "Data size not found", file=sys.stderr)
            #
            # Try to pull the read length
            #
            for rattr in run.xpath("RUN_ATTRIBUTES/*"):
                tag = rattr.findtext("TAG")
                if tag == "actual_read_length":
                    rdata['read_length'] = int(rattr.findtext("VALUE"))
                elif tag == "run":
                    rdata['run_attribute'] = rattr.findtext("VALUE")
            stats = run.xpath("Statistics")
            if stats:
                stats = stats[0]
                sattrib = stats.attrib

                #
                # nreads might lie. cf SRR6263255
                #
                nreads = 0
                for read in stats:
                    rattr = read.attrib
                    # print( rattr, file=sys.stderr)
                    if rattr.has_key('count') and int(rattr['count']) > 0:
                        nreads += 1
                    print(rdata)
                    if rattr.has_key('average') and not 'read_length' in rdata:
                        rdata['read_length'] = float(rattr['average'])
                if nreads > 0:
                    rdata['n_reads'] = nreads
            

            #
            # Calculate estimated disk size based on our data above.
            #
            # Validate first that total_spots * read_size * read_count = total_bases
            #
            # It is possible for all of this to fail. Wrap the entire computation
            # in a try-block; this means upstream computations will need to cope
            # if possible with incomplete size data.
            #
            try:
                if not 'n_reads' in rdata:
                    if exp['library_layout'] == 'PAIRED':
                        rdata['n_reads'] = 2
                    else:
                        rdata['n_reads'] = 1

                if 'read_length' in rdata:
                    calc_bases = rdata['read_length'] * rdata['n_reads'] * rdata['total_spots']
                    err = abs(calc_bases - rdata['total_bases']) / rdata['total_bases']
                    # print( "calc=%d val=%d %f" % (calc_bases, rdata['total_bases'], err), file=sys.stderr)
                    
                    if err > 0.1:
                        print( "Bad size calculation", file=sys.stderr)
                    else:
                        if 'run_id' in rdata:
                            hlen = len(rdata['run_id']) + 50
                        else:
                            hlen = 40
                            # header is 40 ish unless run_id is set; data size is read_length. Two header/data pairs per read,
                            # total_spots of those per file, n_reads files.
                            rdata['estimated_size'] = (hlen + rdata['read_length']) * 2 * rdata['n_reads'] * rdata['total_spots']
                        
                            
            except Exception as e:
                print(f"Failed to compute size data: {e}", file=sys.stderr)
            # print(rdata)
            exp['runs'].append(rdata)


        # just assume one sample
        sample = safe_read(experiment_package, 'SAMPLE', index=0)
        exp['sample_id'] = sample.attrib['accession']
        exp['sample_title'] = safe_read(sample, 'TITLE/text()', index=0)
        exp['sample_description'] = safe_read(sample, 'DESCRIPTION/text()', index=0)
        exp['sample_organism'] = safe_read(sample, 'SAMPLE_NAME/SCIENTIFIC_NAME/text()', index=0)
        exp['sample_taxon'] = safe_read(sample, 'SAMPLE_NAME/TAXON_ID/text()', index=0)
        sample_attrib_tags = safe_read(sample, 'SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE/TAG/text()')
        sample_attrib_vals = safe_read(sample, 'SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE/VALUE/text()')
        sample_attribs = zip(sample_attrib_tags, sample_attrib_vals)
        for sample_attrib in sample_attribs:
            # prepend sample_ for each key related to the sample (unless it's already there)
            if sample_attrib[0].startswith('sample_'):
                    key = sample_attrib[0]
            else:
                key = 'sample_' + sample_attrib[0]
            exp[key] = sample_attrib[1]
        exp_meta.append(exp)

    # create a new run-centric metadata dict
    run_meta = []

    # for each run id in each experiment, create a new record
    for exp in exp_meta:
        for rdata in exp['runs']:
            run = {}
            for k in rdata:
                run[k] = rdata[k]
            for k in exp:
                if k != 'runs':
                    run[k] = exp[k]
            run_meta.append(run)

    return run_meta

def get_runinfo(run_accession):
    """ take sra run accession (like SRR123456)
    return dictionary with keys like: spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path.....
    """
    runinfo_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&rettype=runinfo&id="+run_accession
    # runinfo_url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term="+run_accession
    r = urllib.urlopen(runinfo_url)
    lines = r.read().split("\n")
    keys   = lines[0].split(",")
    values = lines[1].split(",")
    runinfo = dict(zip(keys, values))
    return runinfo

def ftp_download_single_run(run_accession):
    """ use ftp to download one sra file 
    one advantage of this is that it does not leave large file in $HOME/ncbi/public/sra/
    """
    sra_file_url = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s.sra"%(run_accession[:3], run_accession[:6], run_accession, run_accession)
    with open(run_accession+".sra", 'wb') as OUT:
        response = urllib.urlopen(sra_file_url)
        shutil.copyfileobj(response, OUT)

def fastqDumpExistingSraFile(file_name, splitFiles = False):
    """ assumes fastq-dump is on path """
    command = ["fastq-dump"]
    if splitFiles:
        command.append("--split-files")
    command.append(file_name)
    subprocess.call(command, shell=False)

def get_run_ids(run_meta):
    run_ids = []
    for run in run_meta:
        run_ids.append(run['run_id'])
    return run_ids

def download_fastq_files(fasterq_dump_loc, output_dir, metadata, gzip_output=True, splitFiles=False):

    #
    # determine the list of ids to download
    #
    run_ids = [item['run_id'] for item in metadata]

    #
    # Prefetch ids.
    #

    cmd = ['prefetch'] + run_ids
    try:
        result = retry_subprocess_check_output(cmd, 10, 60);

    except subprocess.CalledProcessError as e:
        print(  "Prefetch failed with code " + str(e.returncode), file=sys.stderr)

    # for each run, download the fastq file
    for item in metadata:
        run_id = item['run_id']

        #
        # If not a pacbio run,  use the new fasterq-dump utility (faster, but doesn't have
        # filtering option which results in bigger files)
        #
        # Pacbio requires fastq-dump, and we do not split files for it.
        #
        if item['platform_name'] == 'PACBIO_SMRT':
            fasterq_dump_cmd = ['fastq-dump', '--outdir', output_dir,  run_id]
        else:
            tmpdir = os.getenv("TMPDIR")
            if tmpdir is None:
                tmpdir = "/tmp"
            fasterq_dump_cmd = [fasterq_dump_loc, '-t', tmpdir, '--outdir', output_dir, '-f', run_id]
            if splitFiles: #optional, not needed, --split-3 is the default (see https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)
                fasterq_dump_cmd.append('--split-files')

        try:
            # print( 'executing \'' + str(fasterq_dump_cmd) + '\'', file=sys.stderr)
            result = retry_subprocess_check_output(fasterq_dump_cmd, 5, 60)

        except subprocess.CalledProcessError as e:
            print( "Download failed with code " + str(e.returncode), file=sys.stderr)
            sys.exit(e.returncode)

        # gzip the output; gzip will skip anything already zipped
        if gzip_output:
            fq_files = glob.glob(output_dir + "/*.fastq")
            if fq_files.count > 0:
                gzip_cmd = ["gzip"] + fq_files;
                # print "gzip: " + str(gzip_cmd)
                subprocess.call(gzip_cmd)

    return

def download_sra_data(fasterq_dump_loc, fastq_output_dir, accession_id, metaonly, compress_files, metadata_file, sra_metadata_file):

    # ===== 1. Get the metadata for each run
    metadata = get_accession_metadata(accession_id, sra_metadata_file)

    # ===== 2. Get the fastq files for each run
    if not metaonly:

        download_fastq_files(fasterq_dump_loc, fastq_output_dir, metadata, gzip_output=compress_files)
        # print( 'Fastq Filenames:', file=sys.stderr)
        # print( (glob.glob(fastq_output_dir+'/*.fastq*')), file=sys.stderr)
        for run in metadata:
            run_id = run.get("run_id",None)
            files = glob.glob(os.path.join(fastq_output_dir,"*"+run_id+"*"))
            run['files']=[os.path.basename(f) for f in files]

    # ===== 3. Pack it into the output JSON
    # print( 'Metadata:', file=sys.stderr)
    # print( json.dumps(metadata, indent=1), file=sys.stderr)
    if metadata_file:
        with open(metadata_file, "w") as fp:
            json.dump(metadata, fp, indent=2)



    # ===== 4. Add everything to the workspace (?)

def retry_subprocess_check_output(cmd, n_retries, retry_sleep):
    attempt = 0
    failed = False
    last_error = None

    while attempt < n_retries:
        attempt += 1
        if failed:
            time.sleep(retry_sleep)

        failed = False
        #
        # We run this with stdout/stderr to temp files so we can examine them.
        #
        
        print( "Attempt %d of %d at running %s" % (attempt, n_retries, cmd), file=sys.stderr)
        err = tempfile.TemporaryFile()
        ret =  subprocess.call(cmd, stderr=err)
        print("ret=", ret)
        if ret == 0:
            return

        err.seek(0)
        edata = err.read()
        if re.search("failed to resolve", edata, re.MULTILINE):
            print("Invalid accession", file=sys.stderr)
            print(edata, file=sys.stderr)
            raise RuntimeError("SRA resolution failure")
        print ( "Attempt %d of %d failed at running %s: %s" % (attempt, n_retries, cmd, ret), file=sys.stderr)
        print(edata, file=sys.stderr)
        last_error = ret
        failed = True

    message = "Failed after %s retries running %s" % (attempt, cmd)
    print(message, file=sys.stderr)
    raise Exception(message)


