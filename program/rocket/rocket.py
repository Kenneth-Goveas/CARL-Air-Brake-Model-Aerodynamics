import numpy as np

class Rocket:

    def __init__ (self, **kwargs):

        self._rck_prm = {
            'rck_len' : kwargs['rck_len'],
            'rck_rad' : kwargs['rck_rad'],
            'rck_rgh' : kwargs['rck_rgh']
        }
        self._rck_cmp = kwargs['rck_cmp']

        return

    def getCentPres (self, **kwargs):

        torq_coef = 0
        forc_coef = 0

        for cmp in self._rck_cmp:

            cmp_cent_pres = cmp.getCentPres(**(kwargs | self._rck_prm))
            cmp_lift_coef = cmp.getLiftCoef(**(kwargs | self._rck_prm))

            torq_coef += cmp_cent_pres * cmp_lift_coef
            forc_coef += cmp_lift_coef

        cent_pres = torq_coef / forc_coef

        return cent_pres

    def getLiftCoef (self, **kwargs):
        lift_coef = 0
        for cmp in self._rck_cmp:
            lift_coef += cmp.getLiftCoef(**(kwargs | self._rck_prm))
        return lift_coef

    def getDragCoef (self, **kwargs):
        drag_coef = 0
        for cmp in self._rck_cmp:
            drag_coef += cmp.getDragCoef(**(kwargs | self._rck_prm))
        return drag_coef
