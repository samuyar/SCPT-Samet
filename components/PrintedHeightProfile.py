from math import *

"""
Idee: Es wird die Risshöhe genommen und überprüft, ob es im Interlayer Bereich liegt, falls das der Fall ist werden die Betonwerte
fcm, fct und Ec, angepasst und weitergegeben.

Wo weitergegeben? in Controller (glaube ich)

Die Werte fcm, fct, Ec sollten hier also eine if else Funktion haben und dann abhängig von Interlayer weitergegeben werden

Was muss hier also alles geschehen?

Als Aller erstes müssen die Interlayer Höhen mit Toleranz ermittelt werden.

Als zweites muss die Aktuelle Risshöhe y von der Klasse Geometrie übernommen werden und mit den Toleranzwerten überprüft werden

Danach stopp, neue Zeile in Excel wo true oder False weitergegeben wird damit man überprüfen kann y werte aktualisiert werden etc.
"""

class PrintedHeightProfile:

    def __init__(self, geometrie):

        self.geometrie = geometrie

        print_height = 15 # Druckhöhe je Lage in [mm]
        max_height = 300 # Bauteilhöhe, der Wert ist nicht wichtig, hauptsache er ist >= der eigentlichen Bauteilhöhe
        toleranz = 1 # Toleranz zum Treffen der Interlayers in [mm]

        # erzeugt Liste aller Interlayerhöhen von 0 bis max_height alle print.height, hier also [0,15,30,45,...]
        # + 1 ist nicht die Toleranz hier sondern dafür, dass die letzte Höhe, falls teilbar durch print_height, auch betrachtet wird
        interlayers = list(range(0, max_height + 1, print_height))

        # überprüft, ob y auf einer der interlayers liegt mit toleranz
        is_on_interlayer = any(abs(self.y() - h) <= toleranz for h in interlayers)


    def y(self):
        return self.geometrie.y


