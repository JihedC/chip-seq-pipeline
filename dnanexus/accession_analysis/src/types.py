import urlparse
import pprint
import json
import copy
import dxpy
import common
import re

from itertools import chain
from cached_property import cached_property


SERVER = ''
KEYPAIR = ''

COMMON_METADATA = {
    'lab': 'encode-processing-pipeline',
    'award': 'U41HG006992'
}


logging.basicConfig()
# logging.getLogger("requests").setLevel(logging.WARNING)
logger = logging.getLogger(__name__)
logger.addHandler(dxpy.DXLogHandler())
logger.propagate = False
logger.setLevel(logging.INFO)
logger.info('Logging from the applet is activated')


def dxf_md5(dx_fh):
    logger.debug(
        "in dxf_md5 with handler %s with name %s"
        % (dx_fh, dx_fh.name))
    if 'md5sum' in dx_fh.get_properties():
        md5sum = dx_fh.get_properties()['md5sum']
    else:
        with tempfile.NamedTemporaryFile(delete=True) as tmpfile:
            dxpy.download_dxfile(dx_fh.get_id(), tmpfile.name)
            md5sum = common.md5(tmpfile.name)
        try:
            set_property(dx_fh, {'md5sum': md5sum})
        except Exception as e:
            logger.warning(
                '%s: skipping adding md5sum property to %s.' % (e, dx_fh.name))
    logger.debug('exiting dxf_md5 with %s' % (md5sum))
    return md5sum


def dxf_content_md5(dx_fh):
    logger.debug(
        "in dxf_content_md5 with handler %s with name %s"
        % (dx_fh, dx_fh.name))
    with tempfile.NamedTemporaryFile(delete=True) as tmpfile:
        dxpy.download_dxfile(dx_fh.get_id(), tmpfile.name)
        from magic import from_file
        compressed_mimetypes = [
            "application/x-compress",
            "application/x-bzip2",
            "application/x-gzip"
            ]
        mime_type = from_file(tmpfile.name, mime=True)
        if mime_type in compressed_mimetypes:
            with tempfile.NamedTemporaryFile(delete=True) as uncompressed_tmpfile:
                out, err = common.run_pipe([
                    'cat %s' % (tmpfile.name),
                    'gzip -d'], uncompressed_tmpfile.name)
                md5sum = common.md5(uncompressed_tmpfile.name)
        else:
            md5sum = common.md5(tmpfile.name)
    logger.debug('exiting dxf_content_md5 with %s' % (md5sum))
    return md5sum


def flat(l):
    result = []
    for el in l:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flat(el))
        else:
            result.append(el)
    return result


class Base(object):
    """Initializes class"""

    def __init__(self, data=None):
        if type(data) is dict:
            self.payload = data
        else:
            self.payload = {}

    def make_get_request(self, url):
        if not SERVER and KEYPAIR:
            logger.error('Missing API endpoint and credentials')
        result = common.encoded_get(SERVER + url, keypair=KEYPAIR)
        return result

    def make_post_request(self):
        pass



class EncodedObject(Base):
    """docstring for EncodeObject"""
    def __init__(self, data=None, **kwargs):
        super().__init__(data)
        self.name = self.__class__.__name__
        if kwargs:
            self._encode_repr = self.get_encoded_obj(kwargs)
            logger.info('Initialized {} with unique identifier'.format(self.name))
        else:
            self._encode_repr = {}

    @property
    def uuid(self):
        return self._encode_repr.get('uuid')

    @property
    def id(self):
        return self._encode_repr.get('@id')

    @property
    def accession(self):
        return self._encode_repr.get('accession')

    @property
    def status(self):
        return self._encode_repr.get('status')

    def search_by_field(self, field, value):
        logger.info('Searching ENCODE for {} with {}={}'.format(
                    self.name, field, value))
        url = '/search/?type=%s&%s=%s' % (self.name, field, value)
        result = self.make_get_request(url)
        return result.get('@graph', [])

    def get_by_id(self, field, value, datastore='elasticsearch'):
        logger.info('Getting from ENCODE for {} with {}={}'.format(
                    self.name, field, value))
        if field in ['uuid', 'accession', 'id']:
            result = self.make_get_request(value)
            return result
        else:
            raise ValueError('field value must be unique identifier in ENCODED')

    def get_encoded_obj(self, **kwargs):
        field, value = kwargs.popitem()
        return self.get_by_id(field, value)

    def encode_object(self):
        return self._encode_repr


class File(EncodedObject):
    """docstring for File"""
    def __init__(self, data=None, **kwargs):
        super().__init__(data, kwargs)

    def post(self):
        self.make_post_request()

    def exists(self):
        if self.uuid or self.id:
            return True
        return False

    @property
    def assembly(self):
        return self._encode_repr.get('assembly', None)

    @property
    def format(self):
        return self._encode_repr.get('file_format', None)

    @property
    def biological_replicate_number(self):
        replicate = self._encode_repr.get('replicate', {})
        if replicate:
            return replicate.get('biological_replicate_number')
        return None

    @property
    def md5sum(self):
        return self._encode_repr.get('md5sum')

    @property
    def content_md5sum(self):
        return self._encode_repr.get('content_md5sum')


class Experiment(EncodedObject):
    """docstring for Experiment"""
    def __init__(self, data=None, **kwargs):
        super().__init__(data, kwargs)

    def post(self):
        self.make_post_request()

    @cached_property
    def files(self):
        files = []
        for file_id in self._encode_repr.get('original_files', []):
            files.append(File(id=file_id))
        return files

    @cached_property
    def fastqs(self):
        statuses = ['released', 'in progress', 'uploaded']
        formats = ['fastq', 'fasta']
        fastqs = []
        for file in self.files:
            if file.format in formats and file.status in statuses:
                fastqs.append(file)
        return fastqs

    def fastq_replicates_by_number(self, number):
        return [file for file in self.fastqs if
                file.biological_replicate_number == number]

    def file_by_md5sum(self, md5sum=None, content_md5sum=None):
        files = []
        if md5sum:
            field = 'md5sum'
            value = md5sum
            files = [file for file in self.files if
                     file.md5sum == md5sum and file.status != 'replaced']
        if content_md5sum:
            field = 'content_md5sum'
            value = content_md5sum
            files = [file for file in self.files if
                     file.content_md5sum == content_md5sum and
                     file.status != 'replaced']
        if files:
            try:
                released = [file for file in files if file.status =='released']
                assert len(released) <= 1, 'More than one released file with {}={} found on portal'.format(
                    field, value)
                return released[0]
            except IndexError:
                return files[0]
        return None



class Analysis(object):
    """docstring for Analysis"""
    def __init__(self, analysis_id, **kwargs):
        super().__init__()
        self.analysis_id = analysis_id
        self.analysis = dxpy.describe(analysis_id)
        self.experiment_accession = self.get_experiment_accession()
        self.stages = Stages(self.analysis, **kwargs)
        self.use_content_md5sum = kwargs.get('use_content_md5sum', False)

        try:
            self.experiment = Experiment(accession=self.experiment_accession)
        except AttributeError:
            continue

    def file_at_encode(self, use_content_md5sum):
        match = self.experiment.file_by_md5sum(md5sum=dxf_md5(dx_fh))
        if not match and use_content_md5sum:
            match = self.experiment.file_by_md5sum(
                content_md5sum=dxf_content_md5(dx_fh))
        return match

    def get_experiment_accession(self):
        m_executableName = \
            re.search('(ENCSR[0-9]{3}[A-Z]{3})', self.analysis_object['executableName'])

        m_name = \
            re.search('(ENCSR[0-9]{3}[A-Z]{3})', self.analysis_object['name'])

        if not (m_executableName or m_name):
            logger.error("No experiment accession in name {} or executableName {}.".format(
                         self.analysis_object['name'], self.analysis_object['executableName']))
            return
        elif (m_executableName and m_name):
            executableName_accession = m_executableName.group(1)
            name_accession = m_name.group(1)
            if executableName_accession == name_accession:
                return executableName_accession
            else:
                logger.error(
                    'Different experiment accessions: name {}, executableName {}.'.format(
                        self.analysis_object['name'], self.analysis_object['executableName']))
                return None
        else:
            m = (m_executableName or m_name)
            experiment_accession = m.group(1)
            logger.debug("get_experiment_accession returning {}".format(experiment_accession))
            return experiment_accession

    def encoded_replicate_number(self):
        # this is a fragile way to infer the rep number.  It depends on the name of
        # the mapping analysis.  But since there is nowhere in the anlysis input
        # to put it, the only alternative is to infer it from the list of fastq
        # accessions in the inputs.
        # Probably the best thing is to add that as an input into the analysis,
        # then it would be recorded in the analysis metadata
        logger.debug("in get_encoded_repn with analysis[name] {}".format(
                     self.analysis['name']))

        m_name = re.search(
            'Map ENCSR[0-9]{3}[A-Z]{3} rep(\d+)', self.analysis['name'])

        if not m_name:
            logger.error("Could not infer ENCODED repn from analysis name: {}".format(
                         self.analysis['name']))
            return
        else:
            logger.debug("in get_encoded_repn and found repn to be {}".format(
                         m_name.group(1)))
            encoded_repn = int(m_name.group(1))
            return encoded_repn

    @property
    def is_unreplicated(self):
        return (
            self.analysis['properties'].get('unreplicated_experiment') in ['True', 'true']
            or self.analysis['properties'].get('simplicate_experiment') in ['True', 'true'])

    @property
    def is_unary_control(self):
        return self.analysis['properties'].get('unary_control') in ['True', 'true']

    def is_accessioned(self, file_metadata):
        if file_metadata.get('encode_object'):
            return True
        return False

    def accession_outputs(self, stages, force_patch, force_upload, dryrun):
        files = []
        for stage_name, outputs in stages.items():
            for file_metadata in outputs['output_files']:
                file_id = outputs['stage_metadata']['output'][file_metadata['name']]
                project = outputs['stage_metadata']['project']
                logger.debug('in accessioned_outputs getting handler for file {} in {}'.format(
                    file_id, project))
                dx_file = dxpy.DXFile(file_id, project)
                accessioned_file = self.analysis.file_at_encode(dx_file, use_content_md5sum)
                if accessioned_file:
                    logger.info("Found dx file {} named {} accessioned at ENCODE as {}".format(
                        file_id, dx.name, accessioned_file.accession))
                if accessioned_file and not force_patch and not force_upload:
                    file_metadata.update({'encode_object': accessioned_file})
                    files.append(accessioned_file)
                    continue
                if accessioned_file:
                    if force_patch:
                        logger.info("File already accessioned, but force_patch so patching new metadata")
                    if force_upload:
                        logger.info("File already accessioned, but force_upload so patching new metadata")
                logger.info('Accessioning dx file {} named {}'.format(file_id, dx.name))
                analysis = Analysis(outputs['stage_metadata']['parentAnalysis']).analysis
                dataset_accession = analysis.experiment_accession
                dx_description = dx.describe()


class Stages(object):
    """docstring for Stages"""
    def __init__(self, analysis, **kwargs):
        self.analysis = analysis
        self.fqcheck = kwargs.get('fqcheck', False)

    @cached_property
    def stages(self):
        return [stage['execution'] for stage in self.analysis.get('stages')]

    def scrubbed_stage(self, stage):
        logger.debug('in scrubbed_stage with stage {}'.format(pprint.pformat(stage, depth=3)))
        return stage['input'].get('scrub')

    def stage_by_name(self, stage_name, begins_with=False):
        if begins_with:
            stage = next(
                stage for stage in self.stages
                if stage['name'].startswith(stage_name))
        else:
            stage = next(
                stage for stage in self.stages
                if stage['name'] == stage_name)
        return stage

    def stage_name_by_pattern(self, pattern):
        stages = self.stages
        if not isinstance(stages, list):
            stages = [stages]
        logger.debug('in stage_name_by_pattern with stages {} and pattern {}'.format(
            [stage.get('name') for stage in stages], pattern))
        return next(
            re.match(pattern, stage['name']).group(0)
            for stage in stages
            if re.match(pattern, stage['name'])) or None

    def get_assembly(self, stage_output_tuple):
        stages, output_key = stage_output_tuple
        logger.debug('stages {}'.format(stages))
        logger.debug('in get_assembly with output_key {} and stages:\n{}'.format(
                     output_key, pprint.pformat([stage for stage in stages or []])))
        if not stages:
            return None
        for stage in stages.itervalues():
            output_files = stage.get('output_files')
            if output_files:
                logger.debug('output_files: {}'.format(
                    pprint.pformat([output_file for output_file in output_files])))
                output_file_metadata = next(
                    output_file.get('metadata')
                    for output_file in output_files
                    if output_file.get('name') == output_key)
                return output_file_metadata.get('assembly')


    def pooled_controls(self, rep):
        # this is not surfaced explicitly so must be inferred
        # General:  get the id's of the files actually used for the specified rep
        # and pooled controls.  If the id is the same as the pooled control id then
        # return true.
        # Specifically:
        # starting with the peaks_analysis, get its stages
        # get "ENCODE Peaks" stage
        # get the job id for "ENCODE Peaks"
        # get the control and experiment file ID's for the specified rep
        # find the child jobs of the "ENCODE Peaks" job
        # find the child job where macs2 was run with the experiment file
        # corresponding to the experiment file for this rep
        # get from that child job the file ID of the control
        # if the contol file ID for this rep from ENCODE Peaks is the same as in
        # macs2 then return False else return True
        # Could double-check the log output of ENCODE Peaks to search for the
        # strings "Using pooled controls for replicate 1.", "Using pooled controls
        # for replicate 2." and "Using pooled controls."
        # But there is no corresponding "Not pooling controls"
        # message, so it's underdertermined.

        logger.debug('in pooled_controls with peaks_analysis {}; rep {}'.format(self.analysis['id'], rep))

        ENCODE_Peaks_stage = self.stage_by_name("ENCODE Peaks")

        ENCODE_Peaks_exp_file = \
            ENCODE_Peaks_stage['input']['rep{}_ta'.format(rep)]

        ENCODE_Peaks_ctl_file = \
            ENCODE_Peaks_stage['input']['ctl{}_ta'.format(rep)]

        child_jobs = dxpy.find_jobs(
            parent_job=ENCODE_Peaks_stage['id'],
            name="MACS2",
            project=ENCODE_Peaks_stage['project'],
            describe=True)

        rep_job = next(
            job for job in child_jobs
            if job['describe']['input']['experiment'] == ENCODE_Peaks_exp_file)


        rep_job_ctl_file = rep_job['describe']['input']['control']

        logger.info("Rep{} input control file {}; actually used {}".format(
            rep, ENCODE_Peaks_ctl_file, rep_job_ctl_file))

        if ENCODE_Peaks_ctl_file == rep_job_ctl_file:
            logger.info('Inferred controls not pooled for rep{}'.format(rep))
            return False
        else:
            logger.info('Inferred pooled controls for rep{}'.format(rep))
            return True

    def stage_metadata(self, stage_pattern):
        logger.debug('in get_stage_metadata with analysis {} and stage_pattern {}'.format(self.analysis.get('name'), stage_pattern))
        # unfortunately, very early runs of the IDR pipeline mispelled
        # a stage.  We still need to go back to those runs to harvest QC or for
        # other reasons, so here we just special-case it then over-ride the
        # error.
        try:
            stage_metadata = \
                next(stage for stage in self.stages
                     if re.match(stage_pattern, stage['name']))
        except StopIteration:
            if stage_pattern == "IDR Pooled Pseudoreplicates":
                tmp_metadata = \
                    self.stage_metadata("IDR Pooled Pseudoeplicates")
                tmp_metadata['name'] = "IDR Pooled Pseudoreplicates"
                return tmp_metadata
            else:
                raise
        except:
            raise
        else:
            return stage_metadata

    def gather_encode_fastqs(self):
        # This encoded_replicate_number is the biological_replicate_number at ENCODEd, which
        # needs to be puzzled out from the mapping analysis name or, better, by
        # inferring the rep number from the fastqs actually imported into the
        # analysis
        encoded_repn = self.analysis.encoded_replicate_number()
        experiment_fastqs = self.analysis.experiment.fastq_replicates_by_number(encoded_repn)
        experiment_fastq_accessions = [file.accession for file in experiment_fastqs]
        logger.info('{}: Found accessioned experiment fastqs with accessions {}'.format(
                    self.analysis.experiment.accession, experiment_fastq_accessions))
        return experiment_fastq_accessions

    def gather_input_fastqs(self):
        stage = self.stage_by_name('Gather inputs')
        logger.debug("input_stage['input'] JSON:")
        logger.debug(pprint.pformat(stage['input']))
        input_fastq_accessions = []
        input_fastq_accessions.extend(stage['input']['reads1'])
        logger.debug('reads1 only input_fastq_accessions {}'.format(
                     input_fastq_accessions))

        if stage['input']['reads2']:
            # Coerce into a list here because in earlier versions of the pipeline
            # code, reads2 was just a string.
            reads2 = stage['input']['reads2']
            logger.debug('found reads2 {}'.format(reads2))
            if type(reads2) is list:
                input_fastq_accessions.extend(reads2)
            else:
                input_fastq_accessions.extend([reads2])

        logger.debug('reads1 and reads2 input_fastq_accessions {}'.format(
                     input_fastq_accessions))
        return input_fastq_accessions

    def mapped_read_length(self, raw_mapping_stage, fastq_files):
        crop_length = raw_mapping_stage['output'].get('crop_length')
        if not crop_length or crop_length == 'native':
            logger.warning(
                'crop_length {}. Inferring mapped_read_length from fastqs'.format(crop_length))
            native_lengths = set([fq.get('read_length') for fq in fastq_files])
            try:
                assert (len(native_lengths) == 1 and
                        all([isinstance(rl, int) for rl in native_lengths])), \
                       ('fastqs with different or non-integer read_lengths: {}'.format(
                        [(fq.get('accession'), fq.get('read_length')) for fq in fastq_files]))
            except AssertionError:
                if self.fqcheck:
                    raise
                else:
                    logger.warning(
                        'fastqs with different or non-integer read_lengths: {} But fqcheck is False so ignoring'.format(
                            [(fq.get('accession'), fq.get('read_length')) for fq in fastq_files]))
            except:
                raise
            mapped_read_length = int(next(l for l in native_lengths))
        else:
            mapped_read_length = int(crop_length)
        return mapped_read_length

    def reference_file(self, input_stage):
        # here we get the actual DNAnexus file that was used as the reference
        # need to remain backwards-compatible with analyses that used output_JSON
        input_stage_output = input_stage['output'].get('output_JSON') or input_stage['output']
        reference_file = dxpy.describe(input_stage_output['reference_tar'])

        # and construct the alias to find the corresponding file at ENCODEd
        reference_alias = "dnanexus:" + reference_file.get('id')

        logger.debug('looking for reference file with alias {}'.format(reference_alias))

        reference = File(id=reference_alias)
        assert reference, "Reference file {} not found on Portal".format(reference_alias)
        logger.debug('found reference file {}'.format(reference.get('accession')))
        return reference

    def new_stage_with_pattern(self, pattern):
        stage = {
            self.stage_name_by_pattern(pattern): {
                'input_files': [],
                'output_files': [],
                'qc': [],
                'stage_metadata': {}
            }
        }
        return stage

    def new_input_files(self, fastqs, reference, repn):
        input_files = [
            {'name': 'rep{}_fastqs'.format(repn),
             'derived_from': None,
             'metadata': None,
             'encode_object': fastqs},
            {'name': 'reference',
             'derived_from': None,
             'metadata': None,
             'encode_object': reference}
        ]
        return input_files

    def new_bam_metadata(self, output_type, reference, mapped_read_length):
        bam_metadata = common.merge_dicts({
            'file_format': 'bam',
            'output_type': output_type,
            'assembly': reference.assembly,
            'mapped_read_length': mapped_read_length
        }, COMMON_METADATA)
        return bam_metadata

    def mapping_stage_helper(self, replicate_number):
        # Sets up variables needed for mapping_stage and raw_mapping_stage
        # Returns set of tuples to be unpacked into variables
        encode_fastq_accessions = self.gather_encode_fastqs()
        input_fastq_accessions = self.gather_input_fastqs()
        if self.fqcheck:
            assert set(flat(encode_fastq_accessions)) <= set(flat(input_fastq_accessions)), \
                '{} rep{}: Accessioned experiment fastqs differ from analysis.'.format(
                    self.analysis.experiment.accession, replicate_number)
            return 
        else:
            logger.warning(
                '--fqcheck is False, '
                'so not checking to see if experiment and mapped fastqs match')
        raw_mapping_stage = self.stage_by_name('Map ENCSR')
        filter_qc_stage = self.stage_by_name('Filter and QC')
        scrubbed = self.scrubbed_stage(filter_qc_stage)
        fastq_files = [File(accession=acc) for acc in input_fastq_accessions]
        mapped_read_length = self.mapped_read_length(raw_mapping_stage, fastq_files)
        reference = self.reference_file(self.stage_by_name('Gather inputs'))
        return (fastq_files, mapped_read_length, reference, scrubbed)


    def raw_mapping_stages(self, replicate_number):
        logger.debug(
            'in get_raw_mapping_stages with mapping analysis {} and rep {}'.format(
                mapping_analysis['id'], replicate_number))
        stage_inputs = self.mapping_stage_helper(replicate_number)
        if not stage_inputs:
            return None

        fastq_files, mapped_read_length, reference, scrubbed = stage_inputs

        # Creating mapping_stages
        map_encsr = self.new_stage_with_pattern('Map ENCSR.*')
        if scrubbed:
            filter_and_qc = self.new_stage_with_pattern('Filter and QC.*')
            filter_and_qc['input_files'] = self.new_input_files(fastq_files,
                                                                reference,
                                                                replicate_number)
            filter_and_qc['output_files'] = [
                {'name': 'scrubbed_unfiltered_bam',
                 'derived_from': ['rep{}_fastqs'.format(replicate_number), 'reference'],
                 'metadata': self.new_bam_metadata('unfiltered_alignments',
                                                   reference,
                                                   mapped_read_length)}
            ]
            filter_and_qc['qc'] = [qc]
            mapping_stages = {**map_encsr, **filter_and_qc}
        else:
            map_encsr['input_files'] = self.new_input_files(fastq_files,
                                                            reference,
                                                            replicate_number)
            map_encsr['output_files'] = [
                {'name': 'mapped_reads',
                 'derived_from': ['rep{}_fastqs'.format(replicate_number), 'reference'],
                 'metadata': self.new_bam_metadata('unfiltered_alignments',
                                                   reference,
                                                   mapped_read_length)}]
            map_encsr['qc'] = [qc]
            mapping_stages = map_encsr

        for name, stage in mapping_stages.items():
            if not name.startswith('_'):
                stage.update({'stage_metadata': self.stage_metadata(name)})

        return mapping_stages


    def mapping_stages(self, replicate_number):
        logger.debug(
            'in get_mapping_stages with mapping analysis {} and rep {}'.format(
                self.analysis['id'], replicate_number))
        if not self.analysis:
            logger.warning(
                'get_mapping_stages got empty mapping_analysis, returning None')
            return None

        stage_inputs = self.mapping_stage_helper(replicate_number)
        if not stage_inputs:
            return
        fastq_files, mapped_read_length, reference, scrubbed = stage_inputs
        # Creating three stages of mapping stages
        map_encsr = self.new_stage_with_pattern('Map ENCSR.*')
        filter_and_qc = self.new_stage_with_pattern('Filter and QC.*')
        filter_and_qc['input_files'] = self.new_input_files(fastq_files,
                                                            reference,
                                                            replicate_number)
        filter_and_qc['output_files'] = [
            {'name': 'scrubbed_filtered_bam' if scrubbed else 'filtered_bam',
             'derived_from': ['rep{}_fastqs'.format(replicate_number), 'reference'],
             'metadata': self.new_bam_metadata('alignments',
                                               reference,
                                               mapped_read_length)}
        ]
        filter_and_qc['qc'] = [qc, dup_qc, pbc_qc, filtered_qc, xcor_qc]
        cross_correlation = self.new_stage_with_pattern('Calculate cross-correlation.*')

        # Merge three stages, works only in Python 3.5+
        mapping_stages = {**map_encsr, **filter_and_qc, **cross_correlation}

        for name, stage in mapping_stages.items():
            if not name.startswith('_'):
                stage.update({'stage_metadata': self.stage_metadata(name)})

        return mapping_stages

    def mapping_stages_from_tas(self, tas, reps):
        mapping_jobs = \
            [dxpy.describe(ta['createdBy']['job'])
             for ta in tas]

        mapping_analyses = \
            [dxpy.describe(mapping_job['analysis'])
             for mapping_job in mapping_jobs if mapping_job]

        mapping_stages = []
        for (i, repn) in enumerate(reps):
            mapping_stage = \
                self.mapping_stages(repn)
            if not mapping_stage:
                logger.error('{}: failed to find mapping stages for rep{}'.format(
                             self.analysis['id'], repn))
                return None
            else:
                mapping_stages.append(mapping_stage)

        return mapping_stages


    @cached_property
    def peak_mapping_stages(self):
        # Find the tagaligns actually used as inputs into the analysis
        # Find the mapping analyses that produced those tagaligns
        # Find the filtered bams from those analyses
        # Build the stage dict and return it
        peaks_analysis = self.analysis
        logger.debug('in peak_mapping_stages: peaks_analysis is {} named {}'.format(
            peaks_analysis.get('id'), peaks_analysis.get('name')))

        if peaks_analysis.is_unreplicated:
            reps = [1]
        else:
            reps = [1, 2]


        peaks_stage = self.stage_by_name("ENCODE Peaks")

        tas = [dxpy.describe(peaks_stage['input']['rep%s_ta' % (n)])
               for n in reps]

        return self.mapping_stages_from_tas(tas, reps)


    @cached_property
    def control_mapping_stages(self):
        # Find the control inputs

        peaks_analysis = self.analysis
        logger.debug('in control_mapping_stages with peaks_analysis {}'.format(
            peaks_analysis['id']))

        peaks_stage = self.stage_by_name("ENCODE Peaks")

        # here reps is always 1,2 because the peaks job rep numbers are 1,2 ...
        # these are not the ENCODEd biological_replicate_numbers, which are
        # only known to the analysis via its name, or by going back to ENCODEd and
        # figuring out where the fastqs came from
        reps = sorted([
            re.match("ctl(\d+)_ta", input_key).group(1)
            for input_key in peaks_stage['input'].keys()
            if re.match("ctl(\d+)_ta", input_key)])

        tas = [dxpy.describe(peaks_stage['input']['ctl%s_ta' % (n)])
               for n in reps]

        return self.mapping_stages_from_tas(tas)

    @cached_property
    def histone_peak_stages(self):

        peaks_analysis = self.analysis
        mapping_stages = self.peak_mapping_stages
        control_stages = self.control_mapping_stages

        logger.debug(
            'in get_histone_peak_stages with peaks_analysis {};'.format(peaks_analysis['id']) +
            'experiment {} and len(mapping_stages) {} len(control_stages) {}'.format(
                self.analysis.experiment_accession, len(mapping_stages), len(control_stages)))

        bams = \
            [(mapping_stage,
              'scrubbed_filtered_bam' if any([self.scrubbed_stage(stage) for stage in [mapping_stage.get(stage_name).get('stage_metadata') for stage_name in mapping_stage.keys()]])
              else
              'filtered_bam') for mapping_stage in mapping_stages]

        ctl_bams = \
            [(control_stage,
              'scrubbed_filtered_bam' if any([self.scrubbed_stage(stage) for stage in [control_stage.get(stage_name).get('stage_metadata') for stage_name in control_stage.keys()]])
              else
              'filtered_bam') for control_stage in control_stages]

        assemblies = \
            [self.get_assembly(bam)
             for bam in bams + ctl_bams
             if self.get_assembly(bam)]
        observed_assemblies = set(assemblies)
        assert len(observed_assemblies) == 1, "Different bam assemblies for rep1,2 and control rep1,2 bams: {}".format(assemblies)
        assembly = observed_assemblies.pop()

        if not ctl_bams:
            rep1_ctl = []
            rep2_ctl = []
        else:
            if self.pooled_controls(rep=1):
                rep1_ctl = ctl_bams
            else:
                rep1_ctl = [ctl_bams[0]]

            if not peaks_analysis.is_unreplicated:
                if self.pooled_controls(rep=2):
                    rep2_ctl = ctl_bams
                else:
                    rep2_ctl = [ctl_bams[1]]

        common_file_metadata = copy.copy(COMMON_METADATA)
        common_file_metadata.update({'assembly': assembly})

        narrowpeak_metadata = common.merge_dicts(
            {'file_format': 'bed',
             'file_format_type': 'narrowPeak',
             'file_format_specifications': ['encode:narrowPeak.as'],
             'output_type': 'peaks'},
            common_file_metadata)

        narrowpeak_bb_metadata = common.merge_dicts({
            'file_format': 'bigBed',
            'file_format_type': 'narrowPeak',
            'file_format_specifications': ['encode:narrowPeak.as'],
            'output_type': 'peaks'},
            common_file_metadata)

        replicated_narrowpeak_metadata = common.merge_dicts(
            {'file_format': 'bed',
             'file_format_type': 'narrowPeak',
             'file_format_specifications': ['encode:narrowPeak.as'],
             'output_type': 'replicated peaks'},
            common_file_metadata)

        replicated_narrowpeak_bb_metadata = common.merge_dicts({
            'file_format': 'bigBed',
            'file_format_type': 'narrowPeak',
            'file_format_specifications': ['encode:narrowPeak.as'],
            'output_type': 'replicated peaks'},
            common_file_metadata)

        stable_narrowpeak_metadata = common.merge_dicts(
            {'file_format': 'bed',
             'file_format_type': 'narrowPeak',
             'file_format_specifications': ['encode:narrowPeak.as'],
             'output_type': 'stable peaks'},
            common_file_metadata)

        stable_narrowpeak_bb_metadata = common.merge_dicts({
            'file_format': 'bigBed',
            'file_format_type': 'narrowPeak',
            'file_format_specifications': ['encode:narrowPeak.as'],
            'output_type': 'stable peaks'},
            common_file_metadata)

        fc_signal_metadata = common.merge_dicts({
            'file_format': 'bigWig',
            'output_type': 'fold change over control'},
            common_file_metadata)

        pvalue_signal_metadata = common.merge_dicts({
            'file_format': 'bigWig',
            'output_type': 'signal p-value'},
            common_file_metadata)

        # This is lame because the repns are hard-wired.
        # Needs to be made more general
        rep1_bam = bams[0]
        if not peaks_analysis.is_unreplicated:
            rep2_bam = bams[1]
        # This is lame because it assumes pooled controls is all the controls
        # Needs to be made to get the controls that were actually pooled
        pooled_ctl_bams = ctl_bams
        peak_stages = {
            # derived_from is by name here, will be patched into the file metadata
            # after all files are accessioned
            # derived_from can also be a tuple of (stages,name) to connect to#
            # files outside of this set of stages
            self.stage_name_by_pattern("ENCODE Peaks"): {
                'output_files': [
                    {'name': 'rep1_narrowpeaks',
                     'derived_from': [rep1_bam] + rep1_ctl,
                     'metadata': narrowpeak_metadata},

                    {'name': 'rep1_narrowpeaks_bb',
                     'derived_from': ['rep1_narrowpeaks'],
                     'metadata': narrowpeak_bb_metadata},

                    {'name': 'rep1_pvalue_signal',
                     'derived_from': [rep1_bam] + rep1_ctl,
                     'metadata': pvalue_signal_metadata},

                    {'name': 'rep1_fc_signal',
                     'derived_from': [rep1_bam] + rep1_ctl,
                     'metadata': fc_signal_metadata}

                ] if peaks_analysis.is_unreplicated else [

                    {'name': 'rep1_narrowpeaks',
                     'derived_from': [rep1_bam] + rep1_ctl,
                     'metadata': narrowpeak_metadata},

                    {'name': 'rep2_narrowpeaks',
                     'derived_from': [rep2_bam] + rep2_ctl,
                     'metadata': narrowpeak_metadata},

                    {'name': 'pooled_narrowpeaks',
                     'derived_from': [rep1_bam, rep2_bam] + pooled_ctl_bams,
                     'metadata': narrowpeak_metadata},

                    {'name': 'rep1_narrowpeaks_bb',
                     'derived_from': ['rep1_narrowpeaks'],
                     'metadata': narrowpeak_bb_metadata},

                    {'name': 'rep2_narrowpeaks_bb',
                     'derived_from': ['rep2_narrowpeaks'],
                     'metadata': narrowpeak_bb_metadata},

                    {'name': 'pooled_narrowpeaks_bb',
                     'derived_from': ['pooled_narrowpeaks'],
                     'metadata': narrowpeak_bb_metadata},

                    {'name': 'rep1_pvalue_signal',
                     'derived_from': [rep1_bam] + rep1_ctl,
                     'metadata': pvalue_signal_metadata},

                    {'name': 'rep2_pvalue_signal',
                     'derived_from': [rep2_bam] + rep2_ctl,
                     'metadata': pvalue_signal_metadata},

                    {'name': 'pooled_pvalue_signal',
                     'derived_from': [rep1_bam, rep2_bam] + pooled_ctl_bams,
                     'metadata': pvalue_signal_metadata},

                    {'name': 'rep1_fc_signal',
                     'derived_from': [rep1_bam] + rep1_ctl,
                     'metadata': fc_signal_metadata},

                    {'name': 'rep2_fc_signal',
                     'derived_from': [rep2_bam] + rep2_ctl,
                     'metadata': fc_signal_metadata},

                    {'name': 'pooled_fc_signal',
                     'derived_from': [rep1_bam, rep2_bam] + pooled_ctl_bams,
                     'metadata': fc_signal_metadata}
                ],

                'qc': [],

                'stage_metadata': {}  # initialized below
            },

            # older pipeline versions called this "Overlap" so need a pattern for
            # backward compatibility
            self.stage_name_by_pattern("(Overlap|Final) narrowpeaks"): {
                'output_files': [

                    {'name': 'overlapping_peaks',
                     'derived_from': ['rep1_narrowpeaks', 'rep2_narrowpeaks',
                                      'pooled_narrowpeaks'],
                     'metadata': replicated_narrowpeak_metadata},

                    {'name': 'overlapping_peaks_bb',
                     'derived_from': ['overlapping_peaks'],
                     'metadata': replicated_narrowpeak_bb_metadata}
                ] if not peaks_analysis.is_unreplicated else [
                    {'name': 'overlapping_peaks',
                     'derived_from': ['rep1_narrowpeaks'],
                     'metadata': stable_narrowpeak_metadata},

                    {'name': 'overlapping_peaks_bb',
                     'derived_from': ['overlapping_peaks'],
                     'metadata': stable_narrowpeak_bb_metadata}
                ],

                'qc': [
                    'npeaks_in', 'npeaks_out', 'npeaks_rejected'
                ],

                'stage_metadata': {}  # initialized below
            }
        }

        for stage_name in peak_stages:
            if not stage_name.startswith('_'):
                peak_stages[stage_name].update(
                    {'stage_metadata': self.stage_metadata(stage_name)})

        return [peak_stages]

    def accessioned_outputs(self, stages, use_content_md5sum, tag_outputs=True):
        files = []
        for stage_name, outputs in stages.items():
            for file_metadata in outputs['output_files']:
                file_id = outputs['stage_metadata']['output'][file_metadata['name']]
                project = outputs['stage_metadata']['project']
                logger.debug('in accessioned_outputs getting handler for file {} in {}'.format(
                    file_id, project))
                dx_file = dxpy.DXFile(file_id, project)
                accessioned_file = self.analysis.file_at_encode(dx_file, use_content_md5sum)
                if accessioned_file:
                    logger.info("Found dx file {} named {} accessioned at ENCODE as {}".format(
                        file_id, dx.name, accessioned_file.accession))

                    if tag_outputs:
                        file_metadata.update({'encode_object': accessioned_file})
                    files.append(accessioned_file)
        return files


def accession_histone_analysis_files(peaks_analysis_id, force_patch,
                                     force_upload, dryrun, fqcheck, skip_control,
                                     pipeline_version, use_content_md5sum):
    analysis = Analysis(peaks_analysis_id, fqchek=fqcheck,
                        use_content_md5sum=use_content_md5sum)
    for stages in analysis.stages.control_mapping_stages:
        logger.info('Retrieving accessioned outputs for control mappings')
        analysis.stages.accessioned_outputs(stages, use_content_md5sum)
    for stages in analysis.stages.peak_mapping_stages:
        logger.info('Retrieving accessioned outputs for experiment mappings')
        analysis.stages.accessioned_outputs(stages, use_content_md5sum)
    for stages in analysis.stages.histone_peak_stages:










