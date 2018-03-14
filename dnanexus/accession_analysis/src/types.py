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

    def stage_by_name(self, stage_name):
        peaks_stage = next(
            stage for stage in self.stages
            if stage['name'] == stage_name)

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
