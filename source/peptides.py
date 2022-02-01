
# To assign peptides to each protein fraction.
def peptides(pept, type_protein, values, door, proteinsResult, protein):
    if type_protein in ["omega-gliadin", "gamma-gliadin", "ATIs", "LTP"]:
       if all(value in protein["Description"].lower() for value in values):
           proteinsResult[type_protein] = proteinsResult.get(type_protein, 0) + 1
           door = "close"
    if type_protein == "alpha-gliadin":
        if all(value in protein["Description"].lower() for value in values[0]) or all(value in protein["Description"] for value in values[1]):
            proteinsResult[type_protein] = proteinsResult.get(type_protein, 0) + 1
            door = "close"
    if type_protein in ["HMW", "LMW"]:
        if any(value in protein["Description"].lower() for value in values):
            proteinsResult[type_protein] = proteinsResult.get(type_protein, 0) + 1
            door = "close"
    if type_protein in ["Globulin", "Triticin", "Serpin", "Avenin", "Hordein", "Secalin"]:
        if values in protein["Description"].lower():
            proteinsResult[type_protein] = proteinsResult.get(type_protein, 0) + 1
            door = "close"
    return(door)