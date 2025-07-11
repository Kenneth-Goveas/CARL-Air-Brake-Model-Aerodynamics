import csv
import numpy as np

class Data:

    def __init__ (self, **kwargs):

        csv.register_dialect(
            'carl', delimiter = ',', lineterminator = '\n', escapechar = '\\', quotechar = '"',
            doublequote = False, quoting = csv.QUOTE_NONE
        )

        self._dat_file = open(kwargs['dat_path'], mode = 'w', newline = '')
        self._dat_writ = csv.DictWriter(
                             self._dat_file,
                             fieldnames = [
                                 'Angle of attack (rad)',
                                 'Mach number',
                                 'Brake deployment',
                                 'Center of pressure (m)',
                                 'Lift coefficient',
                                 'Drag coefficient'
                             ],
                             dialect = 'carl'
                         )

        self._dat_writ.writeheader()

        return

    def __del__ (self):
        self._dat_file.close()
        return

    def put (self, **kwargs):

        flw_angl = kwargs['flw_angl']
        flw_mach = kwargs['flw_mach']

        act_brk = kwargs['act_brk']

        cent_pres = kwargs['cent_pres']
        lift_coef = kwargs['lift_coef']
        drag_coef = kwargs['drag_coef']

        self._dat_writ.writerow({
            'Angle of attack (rad)'  : '%+.10e' % flw_angl,
            'Mach number'            : '%+.10e' % flw_mach,
            'Brake deployment'       : '%+.10e' % act_brk,
            'Center of pressure (m)' : '%+.10e' % cent_pres,
            'Lift coefficient'       : '%+.10e' % lift_coef,
            'Drag coefficient'       : '%+.10e' % drag_coef
        })

        return
