import numpy as np

class Base:

    def __init__ (self, **kwargs):

        bas_pos = kwargs['bas_pos']
        bas_rad = kwargs['bas_rad']

        self._cmp_pos = bas_pos
        self._cmp_rad = bas_rad
        self._cmp_area = np.pi * bas_rad**2

        return

    def _getPresArea (self):
        pres_area = self._cmp_area
        return pres_area

    def _getPresDragCoef (self, **kwargs):

        flw_angl = kwargs['flw_angl']
        flw_mach = kwargs['flw_mach']

        if flw_angl < 0.3:
            angl_corr = 1 + 10*flw_angl**2 - 22.2222*flw_angl**3
        else:
            angl_corr = 1.04845 + 1.79105*flw_angl - 3.55519*flw_angl**2 + 1.26691*flw_angl**3

        if flw_mach < 1:
            pres_drag_coef = -(0.12 + 0.13*flw_mach**2)
        else:
            pres_drag_coef = -0.25/flw_mach
        pres_drag_coef *= angl_corr

        return pres_drag_coef

    def getCentPres (self, **kwargs):
        cent_pres = self._cmp_pos
        return cent_pres

    def getLiftCoef (self, **kwargs):
        lift_coef = 0
        return lift_coef

    def getDragCoef (self, **kwargs):

        rck_rad = kwargs['rck_rad']

        rck_area = np.pi * rck_rad**2

        pres_area = self._getPresArea()
        pres_drag_coef = self._getPresDragCoef(**kwargs)

        drag_coef = pres_drag_coef * pres_area / rck_area

        return drag_coef

    def getDesc (self):
        desc = {
            'type'     : 'base',
            'position' : self._cmp_pos,
            'radius'   : self._cmp_rad
        }
        return desc
