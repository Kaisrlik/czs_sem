% Detekce výšky tónu různých hudebních nástrojů na bázi DFT
%   Určete frekvence odpovídající všem půltónům v jedné oktávě mezi tóny c a c',
%   tj. pro tóny c - cis - d - dis - e - f - fis - g - gis - a - ais - h - c'.
%   (Pro řešení této úlohy uvažujte temperované ladění, kdy jednotlivé půltóny
%   jsou od sebe vzdáleny rovnoměrně, tj. poměr jejich kmitočtů je roven 12.
%   odmocnině ze dvou, tj. 2^(1/12). Poměr frekvencí tónů vzdálených o jednu
%   oktávu (např. c a c') je roven vždy 2.)
%   Sledujte v MATLABu, krátkodobá a dlouhodobá spektra signálů různých
%   hudebních nástrojů.
%   Vyberte vhodný krátkodobý segment signálu délky 50 ms a zobrazte vypočítané
%   krátkodobé DFT spektrum analyzovaného segmentu a vysvětlete princip detekce
%   výšky tónu.
%   Opakujte ilustrativně pro segmenty délky 5, 15, 30, 100, 200, 5000 ms
%   a vysvětlete, jak je vhodné volit délku segmentu.
%   Uvažujte vliv doplnění nulami při výpočtu DFT.
%   Signály ke zpracování:
%   tony4.zip - data jsou standardní wav-soubory. K načtení do MATLABu použijte
%   funkci readwav. 
%%


