import bxa.sherpa as bxa
import json

# load the data :
load_pha(1, "c104_3657029905.pha",use_errors=True)
load_bkg(1,"4U_1608-522_3657029905_3c50_bkg_cl_bary_fpm.pha.pi",use_errors=True)
load_rmf(1,"nicer-rmf6s-teamonly-array52.rmf")
load_bkg_rmf(1,"nicer-rmf6s-teamonly-array52.rmf")
load_arf(1, "nicer-arf-consim135o-teamonly-array52.arf")
load_bkg_arf(1, "nicer-arf-consim135o-teamonly-array52.arf")

set_xlog()
set_ylog()
set_stat('cstat')
ignore_id(1,None, 0.3)
ignore_id(1,10.0, None)
notice(0.3, 10.0)
set_xsabund('wilm')

