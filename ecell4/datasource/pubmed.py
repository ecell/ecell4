import re
import numbers

from urllib.request import Request, urlopen, HTTPError  # Python3
from xml.dom import minidom


def read_url(url):
    f = urlopen(url)
    content = f.read().decode('utf-8')
    f.close()
    try:
        f = urlopen(url)
        content = f.read().decode('utf-8')
        f.close()
    except IOError:
        #XXX: say something here
        content = None
    return content

def description(entity):
    entity_id = PubMedDataSource.parse_entity(entity)
    if entity_id is not None:
        entry = []
        src = PubMedDataSource(entity)
        entry.append(('PubMed', '{}'.format(entity_id), ' - '))
        entry.append(('Title', src.data['Title']))
        entry.append(('Author(s)', ', '.join(src.data['AuthorList'])))
        entry.append(('Source', src.data['Source']))
        entry.append(('SO', src.data['SO']))
        entry.append(('DOI', 'https://doi.org/{}'.format(src.data['DOI'])))
        entry.append(('URL', src.link(entity_id)))
        return [entry]

    return []

class PubMedDataSource(object):

    def __init__(self, entity=None):
        self.entity = entity

        if entity is not None:
            entity_id = self.parse_entity(entity)
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={}".format(entity_id)
            data = self.parse_esummary(read_url(url))
            assert len(data) == 1
            self.data = data[0]
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={}&rettype=abstract".format(entity_id)
            abstract = self.parse_abstract(read_url(url))
            if len(abstract) == 0:
                self.abstract = ""
            else:
                assert len(abstract) == 1
                self.abstract = abstract[0]
        else:
            self.data = None
            self.abstract = ""

    @classmethod
    def parse_entity(cls, entity):
        # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000015
        collection = 'pubmed'
        idpttrn = r'\d+'
        uri1 = r'https://www.ncbi.nlm.nih.gov/pubmed/(?P<id>{})'.format(idpttrn)
        uri2 = r'http://identifiers.org/pubmed/(?P<id>{})'.format(idpttrn)
        if isinstance(entity, numbers.Integral) and entity >= 0:
            entity = str(entity)
        if isinstance(entity, str):
            if re.match(r'^{}$'.format(idpttrn), entity) is not None:
                return entity
            mobj = re.match(uri1, entity)
            if mobj is not None:
                return mobj.group('id')
            mobj = re.match(uri2, entity)
            if mobj is not None:
                return mobj.group('id')
        else:
            import ecell4_base
            if isinstance(entity, ecell4_base.core.Species) and entity.has_attribute(collection):
                return cls.parse_entity(entity.get_attribute(collection))
        return None  #XXX: Error

    @classmethod
    def link(cls, entity):
        entity_id = cls.parse_entity(entity)
        assert entity_id is not None
        return 'https://www.ncbi.nlm.nih.gov/pubmed/{}'.format(entity_id)

    @classmethod
    def parse_esummary(cls, esummary):
        retval = []
        doc = minidom.parseString(esummary)
        for entry_node in doc.getElementsByTagName('DocSum'):
            entry = {}
            entry['ID'] = entry_node.getElementsByTagName('Id')[0].firstChild.data
            for item in entry_node.getElementsByTagName('Item'):
                name = item.getAttribute('Name')
                if name in ('Title', 'Volume', 'Issue', 'Pages', 'Source', 'PubDate', 'SO', 'DOI', 'FullJournalName'):
                    if item.firstChild is not None:
                        entry[name] = item.firstChild.data
                elif name == 'AuthorList':
                    entry['AuthorList'] = [author.firstChild.data for author in item.getElementsByTagName('Item') if author.getAttribute('Name') == 'Author']
            retval.append(entry)
        return retval

    @classmethod
    def parse_abstract(cls, efetch):
        retval = []
        doc = minidom.parseString(efetch)
        for node in doc.getElementsByTagName('AbstractText'):
            retval.append(node.firstChild.data)
        return retval

class Formatter(object):

    def __init__(self, entity):
        entity_id = PubMedDataSource.parse_entity(entity)
        if entity_id is None:
            self.src = None
        else:
            self.src = PubMedDataSource(entity)

    @property
    def abstract(self):
        return self.src.abstract

    def __str__(self):
        if self.src is None:
            return None
        authors = ', '.join(self.src.data['AuthorList'])
        year = self.src.data['PubDate'].strip().split(' ')[0]
        data = {key: self.src.data.get(key, '') for key in ('Title', 'FullJournalName', 'Issue', 'Volume', 'Pages', 'DOI', 'ID')}
        return "{Authors}, {Title} {FullJournalName}, {Issue}({Volume}), {Pages}, {Year}. {DOI}. PubMed MPID: {ID}.".format(Authors=authors, Year=year, **data)

    def _ipython_display_(self):
        if self.src is None:
            return
        from IPython.display import display, Markdown
        authors = ', '.join(self.src.data['AuthorList'])
        year = self.src.data['PubDate'].strip().split(' ')[0]
        doi_url = 'https://doi.org/{}'.format(self.src.data['DOI'])
        url = self.src.link(self.src.data['ID'])
        data = {key: self.src.data.get(key, '') for key in ('Title', 'FullJournalName', 'Issue', 'Volume', 'Pages', 'DOI', 'ID')}
        if data["Issue"] != "":
            data["Issue"] = "**{}**".format(data["Issue"])
        if data["FullJournalName"] != "":
            data["FullJournalName"] = "*{}*".format(data["FullJournalName"])
        text = "{Authors}, {Title} {FullJournalName}, {Issue}({Volume}), {Pages}, {Year}. [{DOI}]({DOI_URL}). PubMed PMID: [{ID}]({URL}).".format(Authors=authors, Year=year, DOI_URL=doi_url, URL=url, **data)
        display(Markdown(text))

def citation(entity, formatter=Formatter):
    return formatter(entity)


if __name__ == "__main__":
    print(description("8752322"))
