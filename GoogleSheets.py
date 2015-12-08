__all__ = ['GoogleSheets']

from astropy.table import Table

def GoogleSheets(key, gid, **kwargs):
    url = 'https://docs.google.com/spreadsheets/d/{0}/export?format=csv&gid={1}'.format(key, gid)
    return Table.read(url, format='ascii.csv', **kwargs)

