import urlparse
import pprint
import json
import copy
import dxpy
import common
import re

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
            logger.info('Initialed {} with unique identifier'.format(self.name))
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
        logger.info('Searching ENCODE for {} with {}={}'.format(self.name, field, value))
        url = '/search/?type=%s&%s=%s' % (self.name, field, value)
        result = self.make_get_request(url)
        return result.get('@graph', [])

    def get_by_id(self, field, value, datastore='elasticsearch'):
        logger.info('Getting from ENCODE for {} with {}={}'.format(self.name, field, value))
        if field in ['uuid', 'accession', '@id']
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


class Experiment(EncodedObject):
    """docstring for Experiment"""
    def __init__(self, data=None, **kwargs):
        super().__init__(data, kwargs)
        
    def post(self):
        self.make_post_request()

class Stages(object):
    """docstring for Stages"""
    def __init__(self, arg):
        super(Stages, self).__init__()
        self.arg = arg
        


class Analysis(object):
    """docstring for Analysis"""
    def __init__(self, analysis_id):
        super(Analysis, self).__init__()
        self.analysis_id = analysis_id
        self.analysis = dxpy.describe(analysis_id)
        self.experiment_accession = self.get_experiment_accession()
        try:
            self.experiment = Experiment(accession=self.experiment_accession)
        except AttributeError:
            continue


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
    def is_unreplicated_analysis(self):
        return (
            self.analysis['properties'].get('unreplicated_experiment') in ['True', 'true']
            or self.analysis['properties'].get('simplicate_experiment') in ['True', 'true'])

    @property
    def is_unary_control(self):
        return analysis['properties'].get('unary_control') in ['True', 'true']







