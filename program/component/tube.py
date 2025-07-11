import numpy as np
import scipy as sp

class Tube:

    def __init__ (self, **kwargs):

        tub_pos = kwargs['tub_pos']
        tub_len = kwargs['tub_len']
        tub_rad = kwargs['tub_rad']

        prf_div = 10000
        prf_pos = np.linspace(tub_pos, tub_pos + tub_len, prf_div + 1)
        prf_rad = np.full(prf_div + 1, tub_rad)
        prf_slp = np.full(prf_div + 1, 0)

        self._cmp_pos = prf_pos[0]
        self._cmp_len = prf_pos[prf_div] - prf_pos[0]
        self._cmp_for_rad = prf_rad[0]
        self._cmp_aft_rad = prf_rad[prf_div]
        self._cmp_dif_rad = prf_rad[prf_div] - prf_rad[0]
        self._cmp_for_area = np.pi * prf_rad[0]**2
        self._cmp_aft_area = np.pi * prf_rad[prf_div]**2
        self._cmp_dif_area = np.pi * (prf_rad[prf_div]**2 - prf_rad[0]**2)
        self._cmp_pln_area = 2 * sp.integrate.trapezoid(prf_rad, prf_pos)
        self._cmp_sur_area = 2 * np.pi * sp.integrate.trapezoid(prf_rad * np.sqrt(1 + prf_slp**2), prf_pos)
        self._cmp_vol = np.pi * sp.integrate.trapezoid(prf_rad**2, prf_pos)
        self._cmp_cen = sp.integrate.trapezoid(prf_pos * prf_rad, prf_pos) / sp.integrate.trapezoid(prf_rad, prf_pos) - prf_pos[0]
        self._cmp_tng_inc = np.arctan(prf_slp[prf_div])
        self._cmp_sec_inc = np.arctan((prf_rad[prf_div] - prf_rad[0]) / (prf_pos[prf_div] - prf_pos[0]))

        return

    def _getGljsArea (self):
        gljs_area = self._cmp_pln_area
        return gljs_area

    def _getGljsCentPres (self, **kwargs):
        gljs_cent_pres = self._cmp_pos + self._cmp_cen
        return gljs_cent_pres

    def _getGljsLiftCoef (self, **kwargs):

        flw_angl = kwargs['flw_angl']

        if flw_angl != 0:
            gljs_lift_coef = 1.1 * np.sin(flw_angl)**2 / flw_angl
        else:
            gljs_lift_coef = 0

        return gljs_lift_coef

    def _getFricArea (self):
        fric_area = self._cmp_sur_area
        return fric_area

    def _getFricDragCoef (self, **kwargs):

        flw_angl = kwargs['flw_angl']
        flw_mach = kwargs['flw_mach']

        gas_const = kwargs['gas_const']
        gas_molar = kwargs['gas_molar']
        gas_gamma = kwargs['gas_gamma']
        gas_visc = kwargs['gas_visc']
        gas_temp = kwargs['gas_temp']

        rck_len = kwargs['rck_len']
        rck_rad = kwargs['rck_rad']
        rck_rgh = kwargs['rck_rgh']

        flw_reyn = (flw_mach * rck_len / gas_visc) * np.sqrt(gas_gamma * gas_const * gas_temp / gas_molar)

        if flw_mach < 0.9:
            smt_cmpr_corr = 1 - 0.1*flw_mach**2
        elif flw_mach > 1.1:
            smt_cmpr_corr = (1 + 0.15*flw_mach**2)**(-0.58)
        else:
            smt_cmpr_corr = 5 * ((1.1 - flw_mach) * (1 - 0.1*flw_mach**2) + (flw_mach - 0.9) * (1 + 0.15*flw_mach**2)**(-0.58))

        if flw_mach < 0.9:
            rgh_cmpr_corr = 1 - 0.1*flw_mach**2
        elif flw_mach > 1.1:
            rgh_cmpr_corr = (1 + 0.18*flw_mach**2)**(-1)
        else:
            rgh_cmpr_corr = 5 * ((1.1 - flw_mach) * (1 - 0.1*flw_mach**2) + (flw_mach - 0.9) * (1 + 0.18*flw_mach**2)**(-1))

        if flw_angl < 0.3:
            angl_corr = 1 + 10*flw_angl**2 - 22.2222*flw_angl**3
        else:
            angl_corr = 1.04845 + 1.79105*flw_angl - 3.55519*flw_angl**2 + 1.26691*flw_angl**3

        geom_corr = 1 + rck_rad / rck_len

        if flw_reyn < 10000:
            smt_fric_drag_coef = (1.5*np.log(10000) - 5.6)**(-2)
        else:
            smt_fric_drag_coef = (1.5*np.log(flw_reyn) - 5.6)**(-2)
        smt_fric_drag_coef *= smt_cmpr_corr * angl_corr * geom_corr

        rgh_fric_drag_coef = 0.032 * (rck_rgh / rck_len)**0.2
        rgh_fric_drag_coef *= rgh_cmpr_corr * angl_corr * geom_corr

        fric_drag_coef = max(smt_fric_drag_coef, rgh_fric_drag_coef)

        return fric_drag_coef

    def getCentPres (self, **kwargs):
        cent_pres = self._getGljsCentPres(**kwargs)
        return cent_pres

    def getLiftCoef (self, **kwargs):

        rck_rad = kwargs['rck_rad']

        rck_area = np.pi * rck_rad**2

        gljs_area = self._getGljsArea()
        gljs_lift_coef = self._getGljsLiftCoef(**kwargs)

        lift_coef = gljs_lift_coef * gljs_area / rck_area

        return lift_coef

    def getDragCoef (self, **kwargs):

        rck_rad = kwargs['rck_rad']

        rck_area = np.pi * rck_rad**2

        fric_area = self._getFricArea()
        fric_drag_coef = self._getFricDragCoef(**kwargs)

        drag_coef = fric_drag_coef * fric_area / rck_area

        return drag_coef

    def getDesc (self):
        desc = {
            'type'        : 'tube',
            'position'    : self._cmp_pos,
            'length'      : self._cmp_len,
            'fore radius' : self._cmp_for_rad,
            'aft radius'  : self._cmp_aft_rad
        }
        return desc
