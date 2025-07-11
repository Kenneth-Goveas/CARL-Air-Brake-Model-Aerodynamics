import numpy as np

class Brakes:

    def __init__ (self, **kwargs):

        brk_cnt = kwargs['brk_cnt']
        brk_pos = kwargs['brk_pos']
        brk_area = kwargs['brk_area']

        self._cmp_cnt = brk_cnt
        self._cmp_pos = brk_pos
        self._cmp_area = brk_area

        return

    def _getPresArea (self):
        pres_area = self._cmp_area
        return pres_area

    def _getPresDragCoef (self, **kwargs):

        flw_angl = kwargs['flw_angl']
        flw_mach = kwargs['flw_mach']

        act_brk = kwargs['act_brk']

        if flw_angl < 0.3:
            angl_corr = 1 + 10*flw_angl**2 - 22.2222*flw_angl**3
        else:
            angl_corr = 1.04845 + 1.79105*flw_angl - 3.55519*flw_angl**2 + 1.26691*flw_angl**3

        cnt_corr = self._cmp_cnt

        geom_corr = act_brk

        if flw_mach < 1:
            for_pres_drag_coef = 0.85 + 0.2125*flw_mach**2 + 0.02125*flw_mach**4
        else:
            for_pres_drag_coef = 1.5589 - 0.646/flw_mach**2 + 0.1411/flw_mach**4 + 0.02975/flw_mach**6
        for_pres_drag_coef *= angl_corr * cnt_corr * geom_corr

        if flw_mach < 1:
            aft_pres_drag_coef = 0.12 + 0.13*flw_mach**2
        else:
            aft_pres_drag_coef = 0.25/flw_mach
        aft_pres_drag_coef *= angl_corr * cnt_corr * geom_corr

        pres_drag_coef = for_pres_drag_coef + aft_pres_drag_coef

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
            'type'     : 'brakes',
            'count'    : self._cmp_cnt,
            'position' : self._cmp_pos,
            'area'     : self._cmp_area
        }
        return desc
