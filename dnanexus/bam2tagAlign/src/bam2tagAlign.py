#!/usr/bin/env python
# bam2tagAlign 0.0.1
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# See https://wiki.dnanexus.com/Developer-Portal for documentation and
# tutorials on how to modify this file.
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/

import os, subprocess, shlex
import dxpy

def run_pipe(steps, outfile=None):
    #break this out into a recursive function
    #TODO:  capture stderr
    from subprocess import Popen, PIPE
    p = None
    p_next = None
    first_step_n = 1
    last_step_n = len(steps)
    for n,step in enumerate(steps, start=first_step_n):
        print "step %d: %s" %(n,step)
        if n == first_step_n:
            if n == last_step_n and outfile: #one-step pipeline with outfile
                with open(outfile, 'w') as fh:
                    print "one step shlex: %s to file: %s" %(shlex.split(step), outfile)
                    p = Popen(shlex.split(step), stdout=fh)
                break
            print "first step shlex to stdout: %s" %(shlex.split(step))
            p = Popen(shlex.split(step), stdout=PIPE)
            #need to close p.stdout here?
        elif n == last_step_n and outfile: #only treat the last step specially if you're sending stdout to a file
            with open(outfile, 'w') as fh:
                print "last step shlex: %s to file: %s" %(shlex.split(step), outfile)
                p_last = Popen(shlex.split(step), stdin=p.stdout, stdout=fh)
                p.stdout.close()
                p = p_last
        else: #handles intermediate steps and, in the case of a pipe to stdout, the last step
            print "intermediate step %d shlex to stdout: %s" %(n,shlex.split(step))
            p_next = Popen(shlex.split(step), stdin=p.stdout, stdout=PIPE)
            p.stdout.close()
            p = p_next
    out,err = p.communicate()
    return out,err


@dxpy.entry_point('main')
def main(input_bam, paired_end):

    input_bam_file = dxpy.DXFile(input_bam)

    # The following line(s) download your file inputs to the local file system
    # using variable names for the filenames.

    input_bam_filename = input_bam_file.name
    input_bam_basename = input_bam_file.name.rstrip('.bam')
    dxpy.download_dxfile(input_bam_file.get_id(), input_bam_filename)

    intermediate_TA_filename = input_bam_basename + ".tagAlign"
    if paired_end:
        end_infix = 'PE2SE'
    else:
        end_infix = 'SE'
    final_TA_filename = input_bam_basename + '.' + end_infix + '.tagAlign.gz'

    print subprocess.check_output('ls -l', shell=True)

    # ===================
    # Create tagAlign file
    # ===================
    out,err = run_pipe([
        "bamToBed -i %s" %(input_bam_filename),
        r"""awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}'""",
        "tee %s" %(intermediate_TA_filename),
        "gzip -c"],
        outfile=final_TA_filename)
    print subprocess.check_output('ls -l', shell=True)

    # ================
    # Create BEDPE file
    # ================
    if paired_end:
        final_nmsrt_bam_prefix = input_bam_basename + ".filt.nmsrt.nodup"
        final_nmsrt_bam_filename = final_nmsrt_bam_prefix + ".bam"
        subprocess.check_call(shlex.split("samtools sort -n %s %s" %(input_bam_filename, final_nmsrt_bam_prefix)))

        final_BEDPE_filename = input_bam_basename + ".bedpe.gz"
        out,err = run_pipe([
            "bamToBed -bedpe -mate1 -i %s" %(final_nmsrt_bam_filename),
            "gzip -c"],
            outfile=final_BEDPE_filename)

    print subprocess.check_output('ls -l', shell=True)


    # The following line(s) use the Python bindings to upload your file outputs
    # after you have created them on the local file system.  It assumes that you
    # have used the output field name for the filename for each output, but you
    # can change that behavior to suit your needs.

    tagAlign_file = dxpy.upload_local_file(final_TA_filename)
    if paired_end:
        BEDPE_file = dxpy.upload_local_file(final_BEDPE_filename)

    # The following line fills in some basic dummy output and assumes
    # that you have created variables to represent your output with
    # the same name as your output fields.

    output = {}
    output["tagAlign_file"] = dxpy.dxlink(tagAlign_file)
    if paired_end:
        output["BEDPE_file"] = dxpy.dxlink(BEDPE_file)

    return output

dxpy.run()
