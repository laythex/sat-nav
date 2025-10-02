import sp3
from astropy.time import Time

# from datetime import datetime, timedelta
product = sp3.Product.from_file("ESA0OPSFIN_20230820000_01D_05M_ORB.SP3")

# satellite from SP3 id
satellite = product.satellite_with_id(b"G13")
poly = sp3.narrowed_records_to_piecewise_polynomial(records=satellite.records, window=5, degree=10)

# base+grace_time
c = 1363568400
t = Time(c, format='gps', scale='utc')
itrs = poly(t)
