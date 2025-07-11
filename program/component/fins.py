import numpy as np

class Fins:

    def __init__ (self, **kwargs):

        fin_cnt = kwargs['fin_cnt']
        fin_pos = kwargs['fin_pos']
        fin_rad = kwargs['fin_rad']
        fin_wdt = kwargs['fin_wdt']
        fin_thk = kwargs['fin_thk']
        fin_bas = kwargs['fin_bas']
        fin_tip = kwargs['fin_tip']
        fin_mid_swp = kwargs['fin_mid_swp']
        fin_for_edg = kwargs['fin_for_edg']
        fin_aft_edg = kwargs['fin_aft_edg']

        self._cmp_cnt = fin_cnt
        self._cmp_pos = fin_pos
        self._cmp_rad = fin_rad
        self._cmp_wdt = fin_wdt
        self._cmp_thk = fin_thk
        self._cmp_bas = fin_bas
        self._cmp_tip = fin_tip
        self._cmp_mac = (2/3) * (fin_bas + fin_tip - fin_bas * fin_tip / (fin_bas * fin_tip))
        self._cmp_off = (1/6) * (2 * fin_wdt * np.tan(fin_mid_swp) + fin_bas - fin_tip) * (1 + fin_tip / (fin_bas + fin_tip))
        self._cmp_asp = 4 * fin_wdt / (fin_bas + fin_tip)
        self._cmp_sid_area = (1/2) * fin_wdt * (fin_bas + fin_tip)
        self._cmp_frn_area = fin_wdt * fin_thk
        self._cmp_for_swp = np.arctan(np.tan(fin_mid_swp) + (1/2) * (fin_bas - fin_tip) / fin_wdt)
        self._cmp_mid_swp = fin_mid_swp
        self._cmp_for_edg = fin_for_edg
        self._cmp_aft_edg = fin_aft_edg

        return

    def _getBarrArea (self):
        barr_area = self._cmp_sid_area
        return barr_area

    def _getBarrCentPres (self, **kwargs):

        flw_mach = kwargs['flw_mach']

        flw_beta = np.sqrt(abs(flw_mach**2 - 1))

        ipl_val = np.array([
                      0.25,
                      0,
                      (5.19615*self._cmp_asp - 2) / (10.3923*self._cmp_asp - 3),
                      0.3849*self._cmp_asp / (3.4641*self._cmp_asp - 1)**2
                  ])

        ipl_coef = np.matmul(
                       np.linalg.inv([
                           [1, 0.5, 0.25, 0.125],
                           [0, 1  , 1   , 0.75 ],
                           [1, 2  , 4   , 8    ],
                           [0, 1  , 4   , 12   ]
                       ]),
                       ipl_val
                   )

        if flw_mach < 0.5:
            barr_cent_pres = self._cmp_pos + self._cmp_off + 0.25*self._cmp_mac
        elif flw_mach > 2:
            barr_cent_pres = self._cmp_pos + self._cmp_off \
                             + (3 * flw_beta * self._cmp_asp - 2) / (6 * flw_beta * self._cmp_asp - 3) * self._cmp_mac
        else:
            barr_cent_pres = self._cmp_pos + self._cmp_off \
                             + (ipl_coef[0] + ipl_coef[1] * flw_mach + ipl_coef[2] * flw_mach**2 + ipl_coef[3] * flw_mach**3) * self._cmp_mac

        return barr_cent_pres

    def _getBarrLiftCoef (self, **kwargs):

        flw_angl = kwargs['flw_angl']
        flw_mach = kwargs['flw_mach']

        gas_gamma = kwargs['gas_gamma']

        flw_beta = np.sqrt(abs(flw_mach**2 - 1))

        if self._cmp_cnt < 5:
            cnt_corr = 0.5*self._cmp_cnt
        elif self._cmp_cnt == 5:
            cnt_corr = 0.474*self._cmp_cnt
        elif self._cmp_cnt == 6:
            cnt_corr = 0.4565*self._cmp_cnt
        elif self._cmp_cnt == 7:
            cnt_corr = 0.427*self._cmp_cnt
        elif self._cmp_cnt == 8:
            cnt_corr = 0.405*self._cmp_cnt
        elif self._cmp_cnt > 8:
            cnt_corr = 0.375*self._cmp_cnt

        geom_corr = 1 + self._cmp_rad / (self._cmp_rad + self._cmp_wdt)

        aux = np.array([
                  flw_angl,
                  np.cos(self._cmp_mid_swp)**2,
                  np.pi * self._cmp_asp,
                  np.pi * self._cmp_asp**3,
                  np.sqrt(1 + (0.3 * self._cmp_asp / np.cos(self._cmp_mid_swp))**2),
                  3.01511,
                  0.404959 + 2.67769*gas_gamma,
                  53.4053 + 12.1935*gas_gamma + 17.615*gas_gamma**2,
                  -8.22304,
                  -7.88881 - 20.2855*gas_gamma,
                  -951.314 - 252.618*gas_gamma - 248.211*gas_gamma**2
              ])

        ipl_val = np.array([
                      aux[2] / (1 + aux[4]),
                      0.2 * aux[3] / (aux[1] * aux[4] * (1 + aux[4])**2),
                      aux[5] + aux[6] * aux[0] + aux[7] * aux[0]**2,
                      aux[8] + aux[9] * aux[0] + aux[10] * aux[0]**2
                  ])

        ipl_coef = np.matmul(
                       np.linalg.inv([
                           [1, 0.8, 0.64, 0.512],
                           [0, 1  , 1.6 , 1.92 ],
                           [1, 1.2, 1.44, 1.728],
                           [0, 1  , 2.4 , 4.32 ]
                       ]),
                       ipl_val
                   )

        if flw_mach < 0.8:
            aux = np.array([
                      np.pi * self._cmp_asp,
                      0.5 * flw_beta * self._cmp_asp / np.cos(self._cmp_mid_swp)
                  ])
            barr_lift_coef = aux[0] / (1 + np.sqrt(1 + aux[1]**2))
        elif flw_mach > 1.2:
            aux = np.array([
                      flw_angl,
                      (gas_gamma + 1) * flw_mach**4,
                      4*flw_beta**2,
                      4*flw_beta**4,
                      (gas_gamma + 1) * flw_mach**8,
                      (2*gas_gamma**2 - 7*gas_gamma - 5) * flw_mach**6,
                      10 * (gas_gamma + 1) * flw_mach**4,
                      6*flw_beta**7,
                      2/flw_beta
                  ])
            barr_lift_coef = aux[8] + ((aux[1] - aux[2]) / aux[3]) * aux[0] + ((aux[4] + aux[5] + aux[6] + 8) / aux[7]) * aux[0]**2
        else:
            barr_lift_coef = ipl_coef[0] + ipl_coef[1] * flw_mach + ipl_coef[2] * flw_mach**2 + ipl_coef[3] * flw_mach**3
        barr_lift_coef *= cnt_corr * geom_corr

        return barr_lift_coef

    def _getFricArea (self):
        fric_area = self._cmp_sid_area
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

        cnt_corr = 2*self._cmp_cnt

        geom_corr = 1 + 2 * self._cmp_thk / self._cmp_mac

        if flw_reyn < 10000:
            smt_fric_drag_coef = (1.5*np.log(10000) - 5.6)**(-2)
        else:
            smt_fric_drag_coef = (1.5*np.log(flw_reyn) - 5.6)**(-2)
        smt_fric_drag_coef *= smt_cmpr_corr * angl_corr * cnt_corr * geom_corr

        rgh_fric_drag_coef = 0.032 * (rck_rgh / rck_len)**0.2
        rgh_fric_drag_coef *= rgh_cmpr_corr * angl_corr * cnt_corr * geom_corr

        fric_drag_coef = max(smt_fric_drag_coef, rgh_fric_drag_coef)

        return fric_drag_coef

    def _getPresArea (self):
        pres_area = self._cmp_frn_area
        return pres_area

    def _getPresDragCoef (self, **kwargs):

        flw_angl = kwargs['flw_angl']
        flw_mach = kwargs['flw_mach']

        if flw_angl < 0.3:
            angl_corr = 1 + 10*flw_angl**2 - 22.2222*flw_angl**3
        else:
            angl_corr = 1.04845 + 1.79105*flw_angl - 3.55519*flw_angl**2 + 1.26691*flw_angl**3

        cnt_corr = self._cmp_cnt

        for_geom_corr = np.cos(self._cmp_for_swp)**2

        if self._cmp_aft_edg == 'flat':
            aft_geom_corr = 1
        elif self._cmp_aft_edg == 'round':
            aft_geom_corr = 0.5
        elif self._cmp_aft_edg == 'taper':
            aft_geom_corr = 0

        if self._cmp_for_edg == 'flat':
            if flw_mach < 1:
                for_pres_drag_coef = 0.85 + 0.2125*flw_mach**2 + 0.02125*flw_mach**4
            else:
                for_pres_drag_coef = 1.5589 - 0.646/flw_mach**2 + 0.1411/flw_mach**4 + 0.02975/flw_mach**6
        elif self._cmp_for_edg == 'round':
            if flw_mach < 0.9:
                for_pres_drag_coef = (1 - flw_mach**2)**(-0.417375) - 1
            elif flw_mach > 1:
                for_pres_drag_coef = 1.214 - 0.502/flw_mach**2 + 0.1095/flw_mach**4
            else:
                for_pres_drag_coef = 2.6065 - 1.785*flw_mach
        for_pres_drag_coef *= angl_corr * cnt_corr * for_geom_corr

        if flw_mach < 1:
            aft_pres_drag_coef = 0.12 + 0.13*flw_mach**2
        else:
            aft_pres_drag_coef = 0.25/flw_mach
        aft_pres_drag_coef *= angl_corr * cnt_corr * aft_geom_corr

        pres_drag_coef = for_pres_drag_coef + aft_pres_drag_coef

        return pres_drag_coef

    def getCentPres (self, **kwargs):
        cent_pres = self._getBarrCentPres(**kwargs)
        return cent_pres

    def getLiftCoef (self, **kwargs):

        rck_rad = kwargs['rck_rad']

        rck_area = np.pi * rck_rad**2

        barr_area = self._getBarrArea()
        barr_lift_coef = self._getBarrLiftCoef(**kwargs)

        lift_coef = barr_lift_coef * barr_area / rck_area

        return lift_coef

    def getDragCoef (self, **kwargs):

        rck_rad = kwargs['rck_rad']

        rck_area = np.pi * rck_rad**2

        fric_area = self._getFricArea()
        fric_drag_coef = self._getFricDragCoef(**kwargs)

        pres_area = self._getPresArea()
        pres_drag_coef = self._getPresDragCoef(**kwargs)

        drag_coef = (fric_drag_coef * fric_area + pres_drag_coef * pres_area) / rck_area

        return drag_coef

    def getDesc (self):
        desc = {
            'type'        : 'fins',
            'count'       : self._cmp_cnt,
            'position'    : self._cmp_pos,
            'radius'      : self._cmp_rad,
            'width'       : self._cmp_wdt,
            'thickness'   : self._cmp_thk,
            'base chord'  : self._cmp_bas,
            'tip chord'   : self._cmp_tip,
            'sweep angle' : self._cmp_mid_swp,
            'fore edge'   : self._cmp_for_edg,
            'aft edge'    : self._cmp_aft_edg
        }
        return desc
