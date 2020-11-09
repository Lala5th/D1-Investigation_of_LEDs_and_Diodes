# Description
This is the data and code for analysis taken for Imperials third year experiment D1-Investigation of LEDs and diodes.

# Dataset naming convention

Each dataset has a specific name combined of several parts. The general format is:
setup-Temperature[_DAQ].csv,

`setup` is the setup used, which is made up of 2 parts, a number and a letter. 's' is the fist setup that was used and references the specific cable and resistor, while now it has been replaced mostly with setup 'b' due to its modularity. Setup b has total resistance of ~17 Ohms. Below is a list of some setup pairings, indicating the diodes contained:
- s1: Diode is 1N4001
- b1: Diode is 1N5401
- b2: Diode is 1N4148
- b3: Diode is 1N4001

Any datasets after this will bear the ID of the diode on them. The temperature is either RT (room temperature) or LN (Liquid nitrogen).
