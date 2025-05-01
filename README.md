# algorytmy-zaawansowane

## Instrukcja uruchomienia
> Program nie ma warstwy GUI i działa całkowicie w terminalu.

Aby uruchomić program należy:
1. Otworzyć folder z plikami programu w terminalu (np. `cmd` na Windows).
2. Wykonać komendę `.\studnie.exe -i .\examples\ex1.txt -o .\result.txt`.

Plik wyjściowy znajduje się w tym samym folderze co plik wykonywalny (chyba że została podana inna ścieżka, patrz niżej).

### Argumenty wejściowe
Przy uruchamianiu programu wymagane są dwa argumenty do poprawnego działania programu:
- `-i input_path` (gdzie `input-path` to ścieżka do pliku wejściowego) - podanie pliku wejściowego
- `-o output_path` (gdzie `output-path` to ścieżka do pliku wyjściowego) - określenie gdzie wygenerowany zostanie plik wyjściowy

### Dodatkowe flagi
- `.\studnie.exe -f` wypisze format plików wejściowych i wyjściowych
- `.\studnie.exe -h` wyświetli tekst "pomoc"

## format in

ilość zadań

ilość_studni wielokrotność

(x1, y1) # studnie
(x2, y2)
...
(xn, yn)

(x1, y1) # domy
(x2, y2)
...
(xmn, xmn)

[kolejne zadanie]

## format out

ilość zadań

koszt
(s1, d1)
(s1, d2)
...
(sn, dmn)

[kolejne zadanie]