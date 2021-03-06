import textwrap

class AttributeManager(object):
    """A key-value pair store with optional attribute grouping

    Attributes
    ----------
    attributes: list[list]
        The actual attribute name-value pairs with an optional grouping value
    attribute_dict: dict
        A mapping from attribute name to indices into :attr:`attributes` and
        the group assignments.
    group_dict: dict
        A mapping from group identifier to indices into :attr:`attributes`
    group_counter: int
        The number of attribute groups assigned.

    """

    __slots__ = ('attributes', 'attribute_dict', 'group_dict', 'group_counter')

    def __init__(self, attributes=None):
        """

        Parameters
        ----------
        attributes : Iterable[list], optional
            Attribute name-value pairs with an optional grouping value. If omitted,
            the attribute store will be empty.
        """
        self.attributes = []

        # Internal index attributes
        self.attribute_dict = {}
        self.group_dict = {}
        self.group_counter = 1
        if attributes is not None:
            self._from_iterable(attributes)

    def get_next_group_identifier(self):
        """Retrieve the next un-used attribute group identifier
        and increment the internal counter.

        Returns
        -------
        str
        """

        next_value = self.group_counter
        self.group_counter += 1
        return str(next_value)

    #### Add an attribute to the list and update the lookup tables
    def add_attribute(self, key, value, group_identifier=None):
        """Add an attribute to the list and update the lookup tables

        Parameters
        ----------
        key : str
            The name of the attribute to add
        value : object
            The value of the attribute to add
        group_identifier : str, optional
            The attribute group identifier to use, if any. If not provided,
            no group is assumed.
        """
        items = [key, value]
        if group_identifier is not None:
            items.append(group_identifier)
        self.attributes.append(items)
        index = len(self.attributes) - 1

        #### If there is already one of these, add it to the lists in attribute_dict
        if key in self.attribute_dict:
            self.attribute_dict[key]["indexes"].append(index)
            if group_identifier is not None:
                self.attribute_dict[key]["groups"].append(group_identifier)

        #### Otherwise, create the entry in attribute_dict
        else:
            if group_identifier is not None:
                self.attribute_dict[key] = {"indexes": [
                    index], "groups": [group_identifier]}
            else:
                self.attribute_dict[key] = {"indexes": [index], "groups": []}

        #### If there is a group_identifier, then update the group_dict
        if group_identifier is not None:
            #### If this group already has one or more entries, add to it
            if group_identifier in self.group_dict:
                self.group_dict[group_identifier].append(index)
            #### Else create and entry for the group_identifier
            else:
                self.group_dict[group_identifier] = [index]

    def get_attribute(self, key, group_identifier=None):
        """Get the value or values associated with a given
        attribute key.

        Parameters
        ----------
        key : str
            The name of the attribute to retrieve
        group_identifier : str, optional
            The specific group identifier to return from.

        Returns
        -------
        attribute_value: object or list[object]
            Returns single or multiple values for the requested attribute.
        """
        indices_and_groups = self.attribute_dict[key]
        if group_identifier is None:
            indices = indices_and_groups['indexes']
            if len(indices) > 1:
                return [self.attributes[i][1] for i in indices]
            else:
                return self.attributes[indices[0]][1]
        else:
            groups = indices_and_groups['groups']
            i = groups.index(group_identifier)
            indices = indices['indexes']
            idx = indices[i]
            return self.attributes[idx]

    def replace_attribute(self, key, value, group_identifier=None):
        indices_and_groups = self.attribute_dict[key]
        if group_identifier is None:
            indices = indices_and_groups['indexes']
            if len(indices) > 1:
                raise ValueError("Cannot replace the value of an attribute used multiple times")
            else:
                self.attributes[indices[0]][1] = value
        else:
            raise NotImplementedError()

    def get_by_name(self, name):
        '''Search for an attribute by human-readable name.

        Parameters
        ----------
        name: str
            The name to search for.

        Returns
        -------
        object:
            The attribute value if found or :const:`None`.
        '''
        for attr in self:
            if attr[0].split("|")[-1] == name:
                return attr[1]
        return None

    def clear(self):
        """Remove all content from the store.

        """
        self.attributes = []

        # Internal index attributes
        self.attribute_dict = {}
        self.group_dict = {}
        self.group_counter = 1

    def remove_attribute(self, key, group_identifier=None):
        """Remove the value or values associated with a given
        attribute key from the store.

        This rebuilds the entire store, which may be expensive.

        Parameters
        ----------
        key : str
            The name of the attribute to retrieve
        group_identifier : str, optional
            The specific group identifier to return from.

        """
        indices_and_groups = self.attribute_dict[key]
        if group_identifier is None:
            indices = indices_and_groups['indexes']
            if len(indices) > 1:
                indices = sorted(indices, reverse=True)
                for i in indices:
                    self.attributes.remove(i)
            else:
                self.attributes.remove(indices[0])
        else:
            groups = indices_and_groups['groups']
            i = groups.index(group_identifier)
            indices = indices['indexes']
            idx = indices[i]
            self.attributes.remove(idx)
        attributes = self.attributes
        self.clear()
        self._from_iterable(attributes)

    def has_attribute(self, key):
        """Test for the presence of a given attribute

        Parameters
        ----------
        key : str
            The attribute to test for

        Returns
        -------
        bool
        """
        return key in self.attribute_dict

    def __contains__(self, key):
        return self.has_attribute(key)

    def __getitem__(self, key):
        return self.get_attribute(key)

    def __setitem__(self, key, value):
        if self.has_attribute(key):
            self.replace_attribute(key, value)
        else:
            self.add_attribute(key, value)

    def keys(self):
        return self.attribute_dict.keys()

    def __eq__(self, other):
        if other is None:
            return False
        self_keys = self.keys()
        other_keys = other.keys()
        if len(self_keys) != len(other_keys):
            return False
        for key in self_keys:
            try:
                self_value = self.get_attribute(key)
                other_value = other.get_attribute(key)
            except KeyError:
                return False
            if self_value != other_value:
                return False
        return True

    def __ne__(self, other):
        return not self == other

    def __len__(self):
        return len(self.attributes)

    def __iter__(self):
        return iter(self.attributes)

    def _from_iterable(self, attributes):
        mapping = {}
        for attrib in attributes:
            if len(attrib) == 3:
                # this is a grouped attribute
                if attrib[2] in mapping:
                    remap = mapping[attrib[2]]
                else:
                    remap = self.get_next_group_identifier()
                    mapping[attrib[2]] = remap
                self.add_attribute(attrib[0], attrib[1], remap)
            else:
                self.add_attribute(attrib[0], attrib[1])

    def copy(self):
        """Make a deep copy of the object
        """
        return self.__class__(self.attributes)

    def __repr__(self):
        if len(self) == 0:
            return f"{self.__class__.__name__}([])"
        lines = list(map(str, self.attributes))
        template = "{}([\n{}])"
        return template.format(
            self.__class__.__name__,
            textwrap.indent(',\n'.join(lines), ' ' * 2))
