
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
        self.print_height = 15 # Druckhöhe je Lage in [mm]
        self.max_height = 300 # Bauteilhöhe, der Wert ist nicht wichtig, hauptsache er ist >= der eigentlichen Bauteilhöhe
        self.toleranz = 1 # Toleranz zum Treffen der Interlayers in [mm]

    def y(self):
        return self.geometrie.y

    def is_on_interlayer(self):
        interlayers = list(range(0, self.max_height + 1, self.print_height)) # interlayer höhen von 0 bis max_height je Druckhöhe
        return any(abs(self.y() - h) <= self.toleranz for h in interlayers) # überprüft ob y auf einer der interlayer liegt mit toleranz