import warnings

class IndexBase(object):

    @classmethod
    def from_filename(cls, filename, library=None):
        raise NotImplementedError()

    def offset_for(self, record_label) -> int:
        record = self.record_for(record_label)
        return record.offset

    def record_for(self, record_label):
        record = self.search(record_label)
        if isinstance(record, list):
            warnings.warn(
                f"Multiple records found for {record_label}, using the first")
            record = record[0]
        return record

    def search(self, i, **kwargs):
        raise NotImplementedError()

    def add(self, number, offset, name, analyte, attributes=None):
        raise NotImplementedError()

    def commit(self):
        raise NotImplementedError()

    def check_names_unique(self) -> bool:
        '''Checks that all indexed spectra have unique
        ``spectrum name`` parameters.

        Returns
        -------
        bool:
            Whether the spectrum names in the index are unique.
        '''
        seen = set()
        for record in self:
            if record.name in seen:
                return False
            seen.add(record.name)
        return True
