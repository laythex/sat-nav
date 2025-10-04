import numpy as np
import sp3
from astropy.time import Time
import matplotlib.pyplot as plt

c = 299792458

gps_id = 9
sp3_id = 'G' + '0' * int(gps_id < 10) + f'{gps_id}'

t0 = 732801600
ti = t0 + 3600
tf = ti + 5400

# Считываем телеметрию GRACE
# Время    ID НКА    Псевдодальность C/A
grace_distances = np.loadtxt('GPS1B_2023-03-23_C_04.txt', dtype=np.float64,
                             skiprows=196, usecols=[0, 3, 7])
# Берем необходимый промежуток времени и только первый спутник
filtered_entries = np.where((grace_distances[:, 0] >= ti) & (grace_distances[:, 0] < tf) &
                            (grace_distances[:, 1] == gps_id))[0]
grace_distances = grace_distances[filtered_entries]

# Считываем положения GRACE
# Время    Координата X    Координата Y    Координата Z
grace_locations = np.loadtxt('GNV1B_2023-03-23_C_04.txt', dtype=np.float64,
                             skiprows=196, usecols=[0, 3, 4, 5])
# Берем необходимый промежуток времени
filtered_entries = np.where((grace_locations[:, 0] >= ti) & (grace_locations[:, 0] < tf))[0]
grace_locations = grace_locations[filtered_entries]

# Объединяем данные Grace в одну таблицу (данные по дефолту отстортированы по времени)
# Время    ID НКА    Псевдодальность C/A    Координата X    Координата Y    Координата Z
_, filtered_entries, _ = np.intersect1d(grace_locations[:, 0], grace_distances[:, 0], return_indices=True)
grace_locations = grace_locations[filtered_entries]
grace_data = np.concatenate((grace_distances, grace_locations), axis=1)
grace_data = np.delete(grace_data, 3, axis=1)


# Интерполируем положение GPS
product = sp3.Product.from_file("ESA0OPSFIN_20230820000_01D_05M_ORB.SP3")
satellite = product.satellite_with_id(sp3_id.encode())
gps_poly = sp3.narrowed_records_to_piecewise_polynomial(records=satellite.records, window=5, degree=10)

errors = []
# Проходим по всем точкам телеметрии GRACE
for entry in grace_data:
    # Находим положение GRACE
    grace_time = int(entry[0])
    grace_itrs = np.array(entry[3:6])[np.newaxis].T

    # Находим положение GPS
    gps_time = Time(grace_time + 630763200, format='gps', scale='utc')
    gps_itrs = gps_poly(gps_time).cartesian.xyz.value[np.newaxis].T

    # Поворачиваем координаты GPS
    earth_rotation_rate = 7.2921151467e-5
    propagation_time = entry[2] / c
    phi = -earth_rotation_rate * propagation_time
    rotation_matrix = np.array([[np.cos(phi), np.sin(phi), 0],
                                [-np.sin(phi), np.cos(phi), 0],
                                [0, 0, 1]])
    gps_itrs = rotation_matrix @ gps_itrs

    # Находим разницу расстояния между аппаратами и псевдодальностью
    distance = np.linalg.norm(grace_itrs - gps_itrs)
    pseudodistance = entry[2]
    errors.append((grace_time - t0, abs(distance - pseudodistance) * 1e-3))

plt.scatter(*zip(*errors), s=5)

plt.title(sp3_id)
plt.grid()
plt.xlabel('Время, с')
plt.ylabel('Ошибка дальности, км')

plt.show()
