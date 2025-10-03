import numpy as np
import sp3
from astropy.time import Time
from astropy import units
from astropy.coordinates.representation import CartesianRepresentation
import matplotlib.pyplot as plt

gps_id = 9
t0 = 732801600
ti = t0 + 3600
tf = ti + 5400

# Считываем телеметрию GRACE
# Время    ID НКА    Псевдодальность C/A
grace_distances = np.loadtxt('GPS1B_2023-03-23_C_04.txt', dtype=np.float64,
                             skiprows=196, usecols=[0, 3, 7])
# Берем необходимый промежуток времени и только первый спутник
entries_to_remove = np.where((grace_distances[:, 0] < ti) | (grace_distances[:, 0] >= tf) |
                             (grace_distances[:, 1] != gps_id))[0]
grace_distances = np.delete(grace_distances, entries_to_remove, axis=0)

# Считываем положения GRACE
# Время    Координата X    Координата Y    Координата Z
grace_locations = np.loadtxt('GNV1B_2023-03-23_C_04.txt', dtype=np.float64,
                             skiprows=196, usecols=[0, 3, 4, 5])
# Берем необходимый промежуток времени
entries_to_remove = np.where((grace_distances[:, 0] < ti) | (grace_distances[:, 0] >= tf))[0]
grace_locations = np.delete(grace_locations, entries_to_remove, axis=0)

# Интерполируем положение GPS
product = sp3.Product.from_file("ESA0OPSFIN_20230820000_01D_05M_ORB.SP3")
sp3_id = ('G' + '0' * int(gps_id < 10) + f'{gps_id}').encode()
satellite = product.satellite_with_id(sp3_id)
gps_poly = sp3.narrowed_records_to_piecewise_polynomial(records=satellite.records, window=5, degree=10)

errors = []

# Проходим по всем точкам телеметрии GRACE
for entry_distances in grace_distances:
    # Находим положение GRACE
    grace_time = int(entry_distances[0])
    entry_locations = grace_locations[np.where(grace_locations[:, 0] == entry_distances[0])][0]
    grace_itrs = CartesianRepresentation(entry_locations[1] * units.m,
                                         entry_locations[2] * units.m,
                                         entry_locations[3] * units.m)

    # Находим положение GPS
    gps_time = Time(grace_time + 630763200, format='gps', scale='utc')
    gps_itrs = gps_poly(gps_time).cartesian

    # Находим разницу расстояния между аппаратами и псевдодальностью
    distance = (grace_itrs - gps_itrs).norm().value
    pseudodistance = entry_distances[2]
    errors.append((grace_time - t0, abs(distance - pseudodistance) * 1e-3))

plt.scatter(*zip(*errors), s=5)

plt.title('G' + '0' * int(gps_id < 10) + f'{gps_id}')
plt.grid()
plt.xlabel('Время, с')
plt.ylabel('Ошибка дальности, км')

plt.show()
