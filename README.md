# diophante

Un petit script Python qui résout des équations diophantiennes. Un seul fichier, aucune dépendance.

## Ce qu'il peut faire

- **Équations linéaires à n variables** (`6x + 10y + 14z = 22`) : il donne la solution générale sous forme paramétrique
- **Solutions entières ≥ 0** : il les liste
- **Équation de Pell** (`x² − D·y² = 1` ou `−1`) : il trouve la solution fondamentale et génère les suivantes, même avec D = 61

## Utilisation

Il faut Python 3.8 ou plus.

```bash
python diophante.py "6x + 10y + 14z = 22"
python diophante.py "3x + 5y = 100" -n
python diophante.py "x^2 - 61y^2 = 1" -k 3
```

Sans argument, il se lance en mode interactif et demande l'équation. Une ligne vide pour quitter.

| Option | Effet |
| --- | --- |
| `-n`, `--nonneg` | liste aussi les solutions ≥ 0 |
| `-l`, `--limit` | nombre max de solutions listées (20 par défaut) |
| `-k` | nombre de solutions de Pell à afficher (5 par défaut) |

## Exemples

```
$ python diophante.py "6x + 10y + 14z = 22"
Solution générale (paramètres entiers t1, t2) :
  x = -22 + 7t1 + 3t2
  y = t2
  z = 11 - 3t1 - 2t2

$ python diophante.py "x^2 - 61y^2 = 1" -k 1
(x, y) = (1766319049, 226153980)
```

## Comment ça marche

- **Linéaire :** un algorithme d'Euclide appliqué aux colonnes réduit les coefficients à leur pgcd. La colonne du pgcd donne une solution particulière, et les autres colonnes engendrent toutes les solutions du système homogène.
- **Solutions ≥ 0 :** pour deux variables, une formule directe saute de solution en solution. Au-delà, il énumère et se ramène à ce cas.
- **Pell :** il parcourt la fraction continue de √D jusqu'à tomber sur une réduite qui vaut ±1.

## Ce qu'il ne fait pas

Pas d'équations non linéaires générales (elles sont indécidables en toute généralité), et pas de Pell avec un second membre autre que ±1.
