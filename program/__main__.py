import sys
import multiprocessing as mp

from . import config, data

if __name__ == '__main__':

    cfg_path = sys.argv[1]
    dat_path = sys.argv[2]

    cfg = config.Config(cfg_path = cfg_path)
    dat = data.Data(dat_path = dat_path)

    rck = cfg.getRck()

    gas_prm = cfg.getGasPrm()
    flw_prm_lst = cfg.getFlwPrm()
    act_prm_lst = cfg.getActPrm()

    all_prm_lst = []
    for flw_prm in flw_prm_lst:
        for act_prm in act_prm_lst:
            all_prm_lst.append(gas_prm | flw_prm | act_prm)

    def compute (all_prm):
        cent_pres = rck.getCentPres(**all_prm)
        lift_coef = rck.getLiftCoef(**all_prm)
        drag_coef = rck.getDragCoef(**all_prm)
        dat_pnt = all_prm | {'cent_pres': cent_pres, 'lift_coef': lift_coef, 'drag_coef': drag_coef}
        return dat_pnt

    mp.set_start_method('fork')
    pool = mp.Pool()

    dat_pnt_lst = pool.map(compute, all_prm_lst)

    pool.close()
    pool.join()
    pool.terminate()

    for dat_pnt in dat_pnt_lst:
        dat.put(**dat_pnt)
