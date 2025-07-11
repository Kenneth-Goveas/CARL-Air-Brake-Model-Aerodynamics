import json
import numpy as np

from .. import component, rocket

class Config:

    def __init__ (self, **kwargs):

        cfg_file = open(kwargs['cfg_path'], mode = 'r')
        cfg_data = json.load(cfg_file)
        cfg_file.close()

        self._loadGasPrm(cfg_data)
        self._loadFlwPrm(cfg_data)
        self._loadActPrm(cfg_data)
        self._loadRckPrm(cfg_data)

        return

    def _loadGasPrm (self, cfg_data):

        gas_const = cfg_data['atmosphere']['ideal gas constant']
        gas_molar = cfg_data['atmosphere']['molar mass']
        gas_gamma = cfg_data['atmosphere']['heat capacity ratio']
        gas_visc = cfg_data['atmosphere']['kinematic viscosity']
        gas_temp = cfg_data['atmosphere']['temperature']

        self._gas_prm = {
            'gas_const' : gas_const,
            'gas_molar' : gas_molar,
            'gas_gamma' : gas_gamma,
            'gas_visc'  : gas_visc,
            'gas_temp'  : gas_temp
        }

        return

    def _loadFlwPrm (self, cfg_data):

        flw_angl_rng = cfg_data['airflow']['angle of attack range']
        flw_angl_div = cfg_data['airflow']['angle of attack divisions']
        flw_mach_rng = cfg_data['airflow']['mach number range']
        flw_mach_div = cfg_data['airflow']['mach number divisions']

        flw_angl_rng *= np.pi/180

        self._flw_prm = []
        for flw_angl in np.linspace(0, flw_angl_rng, flw_angl_div + 1):
            for flw_mach in np.linspace(0, flw_mach_rng, flw_mach_div + 1):
                self._flw_prm.append({
                    'flw_angl' : flw_angl,
                    'flw_mach' : flw_mach
                })

        return

    def _loadActPrm (self, cfg_data):

        act_brk_div = cfg_data['actuation']['brake deployment divisions']

        self._act_prm = []
        for act_brk in np.linspace(0, 1, act_brk_div + 1):
            self._act_prm.append({
                'act_brk' : act_brk
            })

        return

    def _loadRckPrm (self, cfg_data):

        rck_cmp = []
        for rck_body_data in cfg_data['rocket']['body']:
            if rck_body_data['type'] == 'nose':
                if rck_cmp:
                    raise ValueError('Invalid body component sequence')
                else:
                    rck_cmp.append(component.Nose(
                        nos_len = rck_body_data['length'],
                        nos_rad = rck_body_data['radius'],
                        nos_ogv = rck_body_data['curvature']
                    ))
            elif rck_body_data['type'] == 'shoulder':
                if not rck_cmp:
                    raise ValueError('Invalid body component sequence')
                else:
                    desc = rck_cmp[-1].getDesc()
                    if desc['type'] != 'tube':
                        raise ValueError('Invalid body component sequence')
                    else:
                        rck_cmp.append(component.Shoulder(
                            shl_pos = desc['position'] + desc['length'],
                            shl_len = rck_body_data['length'],
                            shl_for_rad = desc['aft radius'],
                            shl_aft_rad = rck_body_data['radius']
                        ))
            elif rck_body_data['type'] == 'boattail':
                if not rck_cmp:
                    raise ValueError('Invalid body component sequence')
                else:
                    desc = rck_cmp[-1].getDesc()
                    if desc['type'] != 'tube':
                        raise ValueError('Invalid body component sequence')
                    else:
                        rck_cmp.append(component.Boattail(
                            btl_pos = desc['position'] + desc['length'],
                            btl_len = rck_body_data['length'],
                            btl_for_rad = desc['aft radius'],
                            btl_aft_rad = rck_body_data['radius']
                        ))
            elif rck_body_data['type'] == 'tube':
                if not rck_cmp:
                    raise ValueError('Invalid body component sequence')
                else:
                    desc = rck_cmp[-1].getDesc()
                    if desc['type'] == 'tube':
                        raise ValueError('Invalid body component sequence')
                    else:
                        rck_cmp.append(component.Tube(
                            tub_pos = desc['position'] + desc['length'],
                            tub_len = rck_body_data['length'],
                            tub_rad = desc['aft radius']
                        ))

        rck_rad = 0
        for cmp in rck_cmp:
            desc = cmp.getDesc()
            rck_rad = max(rck_rad, desc['fore radius'], desc['aft radius'])

        desc = rck_cmp[-1].getDesc()
        rck_cmp.append(component.Base(
            bas_pos = desc['position'] + desc['length'],
            bas_rad = desc['aft radius']
        ))

        desc = rck_cmp[-1].getDesc()
        rck_len = desc['position']

        rck_rgh = cfg_data['rocket']['surface']['roughness']

        rck_fins_data = cfg_data['rocket']['fins']

        fins_attached = False
        for cmp in rck_cmp:
            desc = cmp.getDesc()
            if desc['type'] != 'tube':
                continue
            if desc['position'] > rck_fins_data['position']:
                continue
            if desc['position'] + desc['length'] < rck_fins_data['position'] + rck_fins_data['base chord']:
                continue
            rck_cmp.append(component.Fins(
                fin_cnt = rck_fins_data['count'],
                fin_pos = rck_fins_data['position'],
                fin_rad = desc['fore radius'],
                fin_wdt = rck_fins_data['width'],
                fin_thk = rck_fins_data['thickness'],
                fin_bas = rck_fins_data['base chord'],
                fin_tip = rck_fins_data['tip chord'],
                fin_mid_swp = rck_fins_data['sweep angle'] * np.pi / 180,
                fin_for_edg = rck_fins_data['fore edge'],
                fin_aft_edg = rck_fins_data['aft edge']
            ))
            fins_attached = True
            break

        if not fins_attached:
            raise ValueError('Invalid fin placement')

        rck_brakes_data = cfg_data['rocket']['brakes']

        brakes_attached = False
        for cmp in rck_cmp:
            desc = cmp.getDesc()
            if desc['type'] != 'tube':
                continue
            if desc['position'] > rck_brakes_data['position']:
                continue
            if desc['position'] + desc['length'] < rck_brakes_data['position']:
                continue
            rck_cmp.append(component.Brakes(
                brk_cnt = rck_brakes_data['count'],
                brk_pos = rck_brakes_data['position'],
                brk_area = rck_brakes_data['area']
            ))
            brakes_attached = True
            break

        if not brakes_attached:
            raise ValueError('Invalid brake placement')

        self._rck = rocket.Rocket(
                        rck_len = rck_len,
                        rck_rad = rck_rad,
                        rck_rgh = rck_rgh,
                        rck_cmp = rck_cmp
                    )

        return

    def getGasPrm (self):
        return self._gas_prm

    def getFlwPrm (self):
        return self._flw_prm

    def getActPrm (self):
        return self._act_prm

    def getRck (self):
        return self._rck
