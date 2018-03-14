import urlparse
import pprint
import json
import copy
import dxpy
import common
import re

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
        if field in ['uuid', 'accession', 'id']
            result = self.make_get_request(value)
            return result
        else:
            raise ValueError('field value must be unique identifier in ENCODED')

    def get_encoded_obj(self, **kwargs):
        field, value = kwargs.popitem()
        return self.get_by_id(field, value)


class File(EncodedObject):
    """docstring for File"""
    def __init__(self, data=None, **kwargs):
        super().__init__(data, kwargs)


    def post(self):
        self.make_post_request()

    @property
    def format(self):
        return self._encode_repr.get('file_format', None)

    @property
    def biological_replicate_number(self):
        replicate = self._encode_repr.get('replicate', {}):
        if replicate:
            return replicate.get('biological_replicate_number')
        return None
            

class Experiment(EncodedObject):
    """docstring for Experiment"""
    def __init__(self, data=None, **kwargs):
        super().__init__(data, kwargs)
        self.files = []
        self.initialize_files()
        
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
        status = ['released', 'in progress', 'uploaded']
        formats = ['fastq', 'fasta']
        fastqs = []
        for file in self.files:
            if file.format in formats and file.status in statuses:
                fastqs.append(file)
        return fastqs

    def fastq_replicates_by_number(self, number):
        return [file for file in self.fastqs if 
                    file.biological_replicate_number == number]
        


class Analysis(object):
    """docstring for Analysis"""
    def __init__(self, analysis_id):
        super().__init__()
        self.analysis_id = analysis_id
        self.analysis = dxpy.describe(analysis_id)
        self.experiment_accession = self.get_experiment_accession()
        try:
            self.experiment = Experiment(accession=self.experiment_accession)
        except AttributeError:
            continue
        self.


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
    
    @property
    def is_unreplicated(self):
        return (
            self.analysis['properties'].get('unreplicated_experiment') in ['True', 'true']
            or self.analysis['properties'].get('simplicate_experiment') in ['True', 'true'])

    @property
    def is_unary_control(self):
        return analysis['properties'].get('unary_control') in ['True', 'true']

class OutputFile(object):
    """docstring for OutputFile"""
    def __init__(self, arg):
        super().__init__()
        self.arg = arg
        

class Stages(object):
    """docstring for Stages"""
    def __init__(self, analysis, **kwargs):
        self.analysis = analysis
        self.fqchek = kwargs.get('fqchek', False)
        
    @cached_property
    def stages(self):
        return [stage['execution'] for stage in self.analysis.get('stages')]

    def scrubbed_stage(self, stage):
        logger.debug('in scrubbed_stage with stage {}'.format(pprint.pformat(stage, depth=3)))
        return stage['input'].get('scrub')


    def stage_by_name(self, stage_name):
        peaks_stage = next(
            stage for stage in self.stages
            if stage['name'] == stage_name)

    def stage_name_by_pattern(self, pattern):
        stages = self.stages
        if not isinstance(stages, list):
            stages = [stages]
        logger.debug('in get_stage_name with stages {} and pattern {}'.format(
            [stage.get('name') for stage in stages], pattern))
        return next(
            re.match(pattern, stage['name']).group(0)
            for stage in stages
            if re.match(pattern, stage['name'])) or None

    def get_assembly(self, stage_output_tuple):
        stages, output_key = stage_output_tuple
        logger.debug('stages %s' % (stages))
        logger.debug(
            'in get_assembly with output_key %s and stages:\n%s'
            % (output_key, pprint.pformat([stage for stage in stages or []])))
        if not stages:
            return None
        for stage in stages.itervalues():
            output_files = stage.get('output_files')
            if output_files:
                logger.debug('output_files: %s' % (pprint.pformat([output_file for output_file in output_files])))
                output_file_metadata = next(
                    output_file.get('metadata')
                    for output_file in output_files
                    if output_file.get('name') == output_key)
                return output_file_metadata.get('assembly')

    def mapping_stages_from_tas(self, tas):
        mapping_jobs = \
            [dxpy.describe(ta['createdBy']['job'])
             for ta in tas]

        mapping_analyses = \
            [dxpy.describe(mapping_job['analysis'])
             for mapping_job in mapping_jobs if mapping_job]

        mapping_stages = []
        for (i, repn) in enumerate(reps):
            mapping_stage = \
                get_mapping_stages(
                    mapping_analyses[i], keypair, server, fqcheck, repn)
            if not mapping_stage:
                logger.error('%s: failed to find mapping stages for rep%d'
                             % (peaks_analysis['id'], repn))
                return None
            else:
                mapping_stages.append(mapping_stage)

        return mapping_stages

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


    @cached_property
    def peak_mapping_stages(self):
        # Find the tagaligns actually used as inputs into the analysis
        # Find the mapping analyses that produced those tagaligns
        # Find the filtered bams from those analyses
        # Build the stage dict and return it
        
        peaks_analysis = self.analysis
        logger.debug(
            'in get_peak_mapping_stages: peaks_analysis is {} named {}'.format(
                peaks_analysis.get('id'), peaks_analysis.get('name')))

        if peaks_analysis.is_unreplicated:
            reps = [1]
        else:
            reps = [1, 2]


        peaks_stage = self.stage_by_name("ENCODE Peaks")

        tas = [dxpy.describe(peaks_stage['input']['rep%s_ta' % (n)])
               for n in reps]

        return mapping_stages_from_tas(tas)


    @cached_property
    def control_mapping_stages(self):
        # Find the control inputs

        peaks_analysis = self.analysis
        logger.debug(
            'in get_control_mapping_stages with peaks_analysis %s'
            % (peaks_analysis['id']))

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

        return mapping_stages_from_tas(tas)

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
            [get_assembly(bam)
             for bam in bams + ctl_bams
             if get_assembly(bam)]
        observed_assemblies = set(assemblies)
        assert len(observed_assemblies) == 1, "Different bam assemblies for rep1,2 and control rep1,2 bams: %s" % (assemblies)
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
            stage_name_by_pattern("ENCODE Peaks"): {
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
            stage_name_by_pattern("(Overlap|Final) narrowpeaks"): {
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
