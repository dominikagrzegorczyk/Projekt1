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

        def xyz2neup(self, x, y, z):
            '''
            Transformacja współrzędnych ortokartezjańskich (x, y, z) na współrzędne topocentryczne (n,e,up)
            w wyniku przesunięcia początku układu współrzędnych do punktu gdzie znajduję się antena odbiornika,
            a następnie rotacji. Parametry rotacji zależne są od szerokosci (phi) oraz długosc (lam) geodezyjnej anteny.


            Parameters
            ----------

           x : float
               Współrzędna X w układzie orto-kartezjańskim.
           y : float
               Współrzędna Y w układzie orto-kartezjańskim.
           z : float
               Współrzędna Z w układzie orto-kartezjańskim.

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
            R = np.array([[-np.sin(lat) * np.cos(lon), -np.sin(lon), np.cos(lat) * np.cos(lon)],
                          [-np.sin(lat) * np.sin(lon), np.cos(lon), np.cos(lat) * np.sin(lon)],
                          [np.cos(lat), 0, np.sin(lat)]])

            dX = [x, y, z]
            dx = R.T @ dX
            n = dx[0]
            e = dx[1]
            up = dx[2]
            return f"{n:.3f}", f"{e:.20f}", f"{up:.3f}"

        def fl_do_2000(self, lat, lon):
            '''
            Transformacja współrzęnych geodezyjnych szerokosci (phi) i długosci (lam) na współrzędne układu związane z systemem 2000.
            Poprzez odwzorowanie czterostefowe Gaussa-Krugera w pasie 3-stopniowym przy południkach osiowych (15°, 18°, 21°, 24°)
            ponumerowanych odpowiednio: 5, 6, 7 i 8 oraz przy współczynniku skali równym 0.999923.

            Parameters
            ----------
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
                lon = 18
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
            xgk = si + (delta_l ** 2 / 2) * N * np.sin(lat) * np.cos(lat) * (
                        1 + (delta_l ** 2 / 12) * np.cos(lat) ** 2 * (5 - t ** 2 + 9 * n2 + 4 * n2 ** 2) + (
                            delta_l ** 4 / 360) * np.cos(lat) ** 4 * (61 - 58 * t ** 2 + 270 * n2 - 330 * n2 * t ** 2))
            ygk = delta_l * N * np.cos(lat) * (
                        1 + (delta_l ** 2 / 6) * np.cos(lat) ** 2 * (1 - t ** 2 + n2) + (delta_l ** 4 / 120) * np.cos(
                    lat) ** 4 * (5 - 18 * t ** 2 + t ** 4 + 14 * n2 - 58 * n2 * t ** 2))
            mo = 0.999923
            xgk_2000 = xgk * mo
            ygk_2000 = ygk * mo + nr * 1000000 + 500000
            return f"{xgk_2000:.3f}", f"{ygk_2000:.3f}"