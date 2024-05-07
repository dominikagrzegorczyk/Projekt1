import sys
from math import *
import numpy as np

o = object()

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        try:
            model = sys.argv[1]
        except IndexError:
            print('Podaj nazwę modelu (wgrs84, grs80, mars)')
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2


    def xyz2flh(self, X, Y, Z, output='dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny.
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim,

        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        output [STR] - optional, defoulf
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        r = sqrt(X ** 2 + Y ** 2)  # promień
        lat_prev = atan(Z / (r * (1 - self.ecc2)))  # pierwsze przybliilizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001 / 206265:
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev) ** 2)
            h = r / cos(lat_prev) - N
            lat = atan((Z / r) * (((1 - self.ecc2 * N / (N + h)) ** (-1))))
        lon = atan(Y / X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat)) ** 2);
        h = r / cos(lat) - N
        if output == "dec_degree":
            return degrees(lat), degrees(lon), h
        elif output == "dms":
            lat = self.deg2dms(degrees(lat))
            lon = self.deg2dms(degrees(lon))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")

    def flh2xyz(self, lat, lon, h):
        '''
        Algorytm odwrotny do algorytmu hirvonena, polega na transformacji
        współrzędnych geodezyjnych długości, szerokości i wysokości elipsoidalnej (phi, lam, h)
        na współrzędne ortokartezjańskie (x, y, z) z dokoładnoscia do 1 mm.


        Parameters
        ----------
        lat : float
           [stopnie dziesiętne] - szerokosc geodezyjna.
        lon : float
            [stopnie dziesiętne] - Długosc geodezyjna.
        h : float
            [metry] - Wysokosc elipsoidalna w metrach.

        Returns
        -------
        x : float
            Współrzędna X w układzie orto-kartezjańskim.
        y : float
            Współrzędna Y w układzie orto-kartezjańskim.
        z : float
            Współrzędna Z w układzie orto-kartezjańskim.
        '''

        lat = radians(lat)  # fi
        lon = radians(lon)  # lam
        RN = self.a / sqrt(1 - self.ecc2 * sin(lat) ** 2)
        q = RN * self.ecc2 * sin(lat)
        x = (RN + h) * cos(lat) * cos(lon)
        y = (RN + h) * cos(lat) * sin(lon)
        z = (RN + h) * sin(lat) - q
        return f"{x:.3f}", f"{y:.3f}", f"{z:.3f}"

    def xyz2neup(self, x, y, z, x_0, y_0, z_0):
        '''
        Transformacja współrzędnych ortokartezjańskich (x, y, z) na współrzędne topocentryczne (n,e,up)
        w wyniku przesunięcia początku układu współrzędnych do punktu gdzie znajduję się antena odbiornika,
        a następnie rotacji. Parametry rotacji zależne są od szerokosci (phi) oraz długosc (lam) geodezyjnej anteny.


        Parameters
        ----------

        x, y, z : float
        Współrzędne geocentryczne satelitów.

        x_0, y_0, z_0 : float
        Współrzędne geocentryczne anteny.


        Returns
        -------
        n : float
            Północ (Northing) w układzie topocentrycznym
        e : float
            Wschód (Easting) w układzie topocentrycznym.
        up: float
            Góra (Up) w układzie topocentrycznym.

        '''

        lat, lon, h = self.xyz2flh(x, y, z)
        lat = radians(lat)
        lon = radians(lon)
        R = np.array([[-sin(lat) * cos(lon), -sin(lon), cos(lat) * cos(lon)],
                          [-sin(lat) * sin(lon), cos(lon), cos(lat) * sin(lon)],
                          [cos(lat), 0, sin(lat)]])
        xyz_t = np.array([[x - x_0],
                             [y - y_0],
                             [z - z_0]])

        dx = R.T @ xyz_t
        n = dx[0]
        n = np.round(n, decimals=3)
        e = dx[1]
        e = np.round(e, decimals=3)
        up = dx[2]
        up = np.round(up, decimals=3)
        return n,e,up

    def fl_do_2000(self, lat, lon):
        '''
        Transformacja współrzęnych geodezyjnych szerokosci (phi) i długosci (lam) na współrzędne układu związane z systemem 2000.
        Poprzez odwzorowanie czterostefowe Gaussa-Krugera w pasie 3-stopniowym przy południkach osiowych (15°, 18°, 21°, 24°)
        ponumerowanych odpowiednio: 5, 6, 7 i 8 oraz przy współczynniku skali równym 0.999923.

        ---------
        Parameters:
        lat : float
            [stopnie] -  szerokość geodezyjna.
        lon : float
            [stopnie] -  długość geodezyjna.


        Returns
        -------
        x : float
            Współrzędna X w układzie związanym z systemem 2000.
        y : float
            Współrzędna Y w układzie związanym z systemem 2000.
        z : float
            Współrzędna Z w układzie związanym z systemem 2000.

        '''
        if lon >= 13.5 and lon < 16.5:
            lon0 = 15
            nr = 5
        elif lon >= 16.5 and lon < 19.5:
            lon0 = 18
            nr = 6
        elif lon >= 19.5 and lon < 22.5:
            lon0 = 21
            nr = 7
        elif lon >= 22.5 and lon <= 25.5:
            lon0 = 24
            nr = 8

        lat = radians(lat)
        lon = radians(lon)
        lon0 = radians(lon0)

        b2 = self.a ** 2 * (1 - self.ecc2)
        e2_poch = (self.a ** 2 - b2) / b2
        delta_l = lon - lon0
        t = tan(lat)
        n2 = e2_poch * np.cos(lat) ** 2
        N = self.a / sqrt(1 - self.ecc2 * sin(lat) ** 2)
        A0 = 1 - self.ecc2 / 4 - 3 * self.ecc2 ** 2 / 64 - 5 * self.ecc2 ** 3 / 256
        A2 = (3 / 8) * (self.ecc2 + self.ecc2 ** 2 / 4 + 15 * self.ecc2 ** 3 / 128)
        A4 = (15 / 256) * (self.ecc2 ** 2 + 3 * self.ecc2 ** 3 / 4)
        A6 = (35 * self.ecc2 ** 3) / 3072
        si = self.a * (A0 * lat - A2 * np.sin(2 * lat) + A4 * np.sin(4 * lat) - A6 * np.sin(6 * lat))
        xgk = si + (delta_l ** 2 / 2) * N * np.sin(lat) * np.cos(lat) * (1 + (delta_l ** 2 / 12) * np.cos(lat) ** 2 * (5 - t ** 2 + 9 * n2 + 4 * n2 ** 2) + (delta_l ** 4 / 360) * np.cos(lat) ** 4 * (61 - 58 * t ** 2 + 270 * n2 - 330 * n2 * t ** 2))
        ygk = delta_l * N * np.cos(lat) * (1 + (delta_l ** 2 / 6) * np.cos(lat) ** 2 * (1 - t ** 2 + n2) + (delta_l ** 4 / 120) * np.cos(lat) ** 4 * (5 - 18 * t ** 2 + t ** 4 + 14 * n2 - 58 * n2 * t ** 2))
        mo = 0.999923
        xgk_2000 = xgk * mo
        ygk_2000 = ygk * mo + nr * 1000000 + 500000
        return f"{xgk_2000:.3f}", f"{ygk_2000:.3f}"

    def fl_do_1992(self, lat, lon, lon0=19):
        '''
        Transformacja współrzęnych geodezyjnych szerokosci (phi) i długosci (lam) na współrzędne układu związane z systemem 1992.
        Poprzez odwzorowanie Gaussa-Krugera w pasie 10-stopniowym przy południku osiowym równum 19°
        oraz przy współczynniku skali na południku osiowym równym 0.9993.


        Parameters
        ----------
       lat : float
           [stopnie] - szerokość geodezyjna.
       lon : float
           [stopnie] - długość geodezyjna.
       lon0 : float
           [stopnie] - południk osiowy (wartosc stała).



        Returns
        -------
        x : float
            Współrzędna X w układzie związanym z systemem 2000.
        y : float
            Współrzędna Y w układzie związanym z systemem 2000.
        z : float
            Współrzędna Z w układzie związanym z systemem 2000.

        '''
        lat = radians(lat)
        lon = radians(lon)
        lon0 = radians(lon0)

        b2 = self.a ** 2 * (1 - self.ecc2)
        e2_poch = (self.a ** 2 - b2) / b2
        delta_l = lon - lon0
        t = tan(lat)
        n2 = e2_poch * np.cos(lat) ** 2
        N = self.a / sqrt(1 - self.ecc2 * sin(lat) ** 2)
        A0 = 1 - self.ecc2 / 4 - 3 * self.ecc2 ** 2 / 64 - 5 * self.ecc2 ** 3 / 256
        A2 = (3 / 8) * (self.ecc2 + self.ecc2 ** 2 / 4 + 15 * self.ecc2 ** 3 / 128)
        A4 = (15 / 256) * (self.ecc2 ** 2 + 3 * self.ecc2 ** 3 / 4)
        A6 = (35 * self.ecc2 ** 3) / 3072
        si = self.a * (A0 * lat - A2 * np.sin(2 * lat) + A4 * np.sin(4 * lat) - A6 * np.sin(6 * lat))
        xgk = si + (delta_l ** 2 / 2) * N * np.sin(lat) * np.cos(lat) * (
                    1 + (delta_l ** 2 / 12) * np.cos(lat) ** 2 * (5 - t ** 2 + 9 * n2 + 4 * n2 ** 2) + (
                        delta_l ** 4 / 360) * np.cos(lat) ** 4 * (61 - 58 * t ** 2 + 270 * n2 - 330 * n2 * t ** 2))
        ygk = delta_l * N * np.cos(lat) * (
                    1 + (delta_l ** 2 / 6) * np.cos(lat) ** 2 * (1 - t ** 2 + n2) + (delta_l ** 4 / 120) * np.cos(
                lat) ** 4 * (5 - 18 * t ** 2 + t ** 4 + 14 * n2 - 58 * n2 * t ** 2))
        mo = 0.9993
        xgk_92 = xgk * mo - 5300000
        ygk_92 = ygk * mo + 500000
        return f"{xgk_92:.3f}", f"{ygk_92:.3f}"


if __name__ == "__main__":
    # utworzenie obiektu
    geo = Transformacje(model="grs80")
    print(sys.argv)
    # dane XYZ geocentryczne
    X = 3664940.50
    Y = 1409153.590
    Z = 5009571.170
    phi, lam, h = geo.xyz2flh(X, Y, Z)
    X_new, Y_new, Z_new = geo.flh2xyz(phi, lam, h)
    x_0 = 1
    y_0 = 1
    z_0 = 1
    n, e, up = geo.xyz2neup(X, Y, Z, x_0, y_0, z_0)
    X_2000, Y_2000 = geo.fl_do_2000(phi, lam)
    X_1992, Y_1992 = geo.fl_do_1992(phi, lam)

    input_file_path = sys.argv[-1]
    if '--header_lines' in sys.argv:
        header_lines = int(sys.argv[3])

    if '--xyz2flh' in sys.argv and '--flh2xyz' in sys.argv and '--xyz2neup' in sys.argv and '--fl_do_2000' in sys.argv and '--fl_do_1992':
        print('mozezz podac tylko jedna flage')
    elif '--xyz2flh' in sys.argv:
        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            lines = lines[header_lines:]
            coords_flh = []
            for line in lines:
                line = line.strip()
                x_str, y_str, z_str = line.split(',')
                x, y, z = (float(x_str), float(y_str), float(z_str))
                f, l, h = geo.xyz2flh(x, y, z)
                coords_flh.append([f, l, h])

        with open('results_xyz2flh.txt', 'w') as f:
            f.write('     phi[deg]           lam[deg]            h[m] \n')
            for coords in coords_flh:
                coords_flh_line = ' '.join([str(coord) for coord in coords])
                f.write(coords_flh_line + '\n')

    elif '--flh2xyz' in sys.argv:
        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            lines = lines[1:]
            coords_xyz_new = []
            for line in lines:
                line = line.strip()
                f_str, l_str, h_str = line.split()
                f, l, h = (float(f_str), float(l_str), float(h_str))
                X_new, Y_new, Z_new = geo.flh2xyz(f, l, h)
                coords_xyz_new.append([X_new, Y_new, Z_new])

        with open('results_flh2xyz.txt', 'w') as f:
            f.write('   X[m]        Y[m]        Z[m]\n')
            for coords in coords_xyz_new:
                coords_xyz_new_line = ' '.join([str(coord) for coord in coords])
                f.write(coords_xyz_new_line + '\n')


    elif '--xyz2neup' in sys.argv:
        with open(input_file_path,'r') as f:
            lines = f.readlines()
            lines = lines[header_lines:]
            # coords_neup = []
            coords_neup = []
            for line in lines:
                line = line.strip()
                x_str, y_str, z_str = line.split(',')
                x, y,z = (float(x_str), float(y_str), float(z_str))
                x_0,y_0,z_0 = [float(coord) for coord in sys.argv[-4:-1]]
                n,e,up = geo.xyz2neup(x,y,z,x_0,y_0,z_0)
                coords_neup.append([n,e,up])


        with open('results_neupxyz.txt', 'w') as f:
            f.write('n[m]     e[m]     up[m] \n')
            for coords in coords_neup:
                coords_neup_line = ' '.join([str(coord) for coord in coords])
                f.write(coords_neup_line + '\n')


    elif '--fl_do_2000' in sys.argv:
        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            lines = lines[1:]
            coords_xy2000 = []
            for line in lines:
                line = line.strip()
                f_str, l_str, h_str = line.split()
                f, l, h = (float(f_str), float(l_str), float(h_str))
                X_2000, Y_2000 = geo.fl_do_2000(f, l)
                coords_xy2000.append([X_2000, Y_2000])

        with open('results_fl_do_2000.txt', 'w') as f:
            f.write('   X[m]        Y[m]\n')
            for coords in coords_xy2000:
                coords_xy2000_line = ' '.join([str(coord) for coord in coords])
                f.write(coords_xy2000_line + '\n')

    elif '--fl_do_1992' in sys.argv:
        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            lines = lines[1:]
            coords_xy1992 = []
            for line in lines:
                line = line.strip()
                f_str, l_str, h_str = line.split()
                f, l, h = (float(f_str), float(l_str), float(h_str))
                X_1992, Y_1992 = geo.fl_do_1992(f, l)
                coords_xy1992.append([X_1992, Y_1992])

        with open('results_fl_do_1992.txt', 'w') as f:
            f.write('   X[m]       Y[m]\n')
            for coords in coords_xy1992:
                coords_xy1992_line = ' '.join([str(coord) for coord in coords])
                f.write(coords_xy1992_line + '\n')
