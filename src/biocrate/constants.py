import periodictable
MAP_SIZE = 1073741824 # 1GB
elem2num = lambda x: periodictable.elements.symbol(x[:1].upper() + x[1:].lower()).number
elem2mass = lambda x: periodictable.elements.symbol(x[:1].upper() + x[1:].lower()).mass
main_enum = {"N": 0, "CA": 1, "C": 2, "O": 3}
