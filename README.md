# Projekt1 Informatyka Geodezyjna 

## Cel i funkcje programu
Program służy do tranformacji współrzędnych pomiędzy poniższymi układami:
- X,Y,Z do fi, lambda, H (xyz2flh)
- fi, lambda, H do X,Y,Z (flh2xyz)
- X,Y,Z na N, E, Up (xyz2neup)
- fi, lambda do układu PL-2000(fldo2000)
- fi, lambda do układu PL-1992(fldo1992)

Program obsługuje elipsoidy:  
**WGS'84**  
**GRS'80**  
**mars**

## Wymogi techniczne
Do poprawnego działania programu należy spełnić warunki:
- Zainstalowany Python w wersji **3.12**
- Zainstalowana bibioteka **numpy** (pip install numpy)
- System operacyjny **Windows 11** mający **CLI**

## Użycie programu

Przed skorzystaniem z programu należy pobrać zawartość repozytorium, a następnie uruchomić interfejs wiersza poleceń (**CLI**) w wcześniej skopiowanym folderze.

Aby program działał należy użyć:

### CLI z flagami 

Umożiwia to wywoływanie programu przy pomocy flag, wszystkie wymagane informacje powinny być podawane w jednej linii.

Do poprawnego uruchomienia programu należy podać:
-flage do określenia modelu elipsoidy **--model** (wgs84,grs80,mars) 
>--model wgs84	--model grs80	--model mars
- flage do okreslenia liczby linijek kodu nagłówka, które użytkownik chce usunąć **--header_lines** 
- flage do określenia funkcji, którą użytkownik chce wykonać
>--xyz2flh	--flh2xyz	--xyz2neup	--fl_do_2000	--fl_do_1992
- plik z którego użytkownik chcę aby zostały odczytane współrzędne początkowe **file.txt**

## Przykładowe wywołanie programu:
1. Przejście do CLI z skopiowanego folderu zawierący program
2. Wywołanie przykładowego programu 
>python projekt1.py --model grs80 --header_lines 4 --xyz2flh wsp_inp.txt 

### Uwaga
Przy wywoływaniu programu z flagą **--xyz2neup** konieczne jest podanie współrzędnych środka układu odniesinia
>python projekt.py --model grs80 --header_lines 4 --xyz2neup 3664940.500 1409153.590 5009571.170 file.txt

Plik z wynikami zostanie zapisany do folderu, z którego został wywołany program. Plik będzie zawierać jeden wiersz nagłówka oraz wartości będą rozdzielone " "(whitespace).

#### Przykład (**results_xyz2flh.txt**)

>phi[deg]      lam[deg]      h[m]  
>52.09727222 21.03153333 141.399  
>52.09727216 21.03153314 141.400  
>52.09727212 21.03153296 141.403  
>52.09727209 21.03153277 141.408  
>52.09727209 21.03153323 141.410  
>52.09727212 21.03153318 141.402  
>52.09727207 21.03153300 141.407  
>52.09727206 21.03153281 141.411  
>52.09727212 21.03153325 141.407  
>52.09727214 21.03153318 141.405  
>52.09727210 21.03153332 141.408  
>52.09727215 21.03153318 141.406  


W pliku z danymi początkowymi X,Y,Z współrzędne powinny być odzielone "**,**" 
>wspolrzednaX,wspolrzednaY,wspolrzednaZ

#### Przykład (**wsp_inp.txt**)

>3664940.500,1409153.590,5009571.170  
>3664940.510,1409153.580,5009571.167  
>3664940.520,1409153.570,5009571.167  
>3664940.530,1409153.560,5009571.168  
>3664940.520,1409153.590,5009571.170  
>3664940.514,1409153.584,5009571.166  
>3664940.525,1409153.575,5009571.166  
>3664940.533,1409153.564,5009571.169  
>3664940.515,1409153.590,5009571.170  
>3664940.514,1409153.584,5009571.169  
>3664940.515,1409153.595,5009571.169  
>3664940.513,1409153.584,5009571.171  
