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
            print('Podaj nazwę modelu')
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
        self.ecc = sqrt(2 * self.flat - self.flat ** 2)  # eccentricity  WGS84:0.0818191910428
        self.ecc2 = (2 * self.flat - self.flat ** 2)  # eccentricity**2


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
        na współrzędne ortokartezjańskie (x, y, z).



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
