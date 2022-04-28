# To assign peptides to each protein fraction.

def peptidesv2(pept, type_protein, values, door, proteinsResult, protein, NCBI_Acc_annotation, listPeptides):
    description = NCBI_Acc_annotation[protein["NCBI_Acc"]]
    if type_protein in ["omega-gliadin", "gamma-gliadin", "ATIs", "LTP"]:
       if all(value in description.lower() for value in values):
           if pept not in listPeptides[type_protein]:
                proteinsResult[type_protein] = proteinsResult.get(type_protein, 0) + 1
                door = "close"
                listPeptides[type_protein].append(pept)
    if type_protein == "alpha-gliadin":
        if all(value in description.lower() for value in values[0]) or all(value in description.lower() for value in values[1]):
            if pept not in listPeptides[type_protein]:
                proteinsResult[type_protein] = proteinsResult.get(type_protein, 0) + 1
                door = "close"
                listPeptides[type_protein].append(pept)
    if type_protein in ["HMW", "LMW"]:
        if any(value in description.lower() for value in values):
            if pept not in listPeptides[type_protein]:
                proteinsResult[type_protein] = proteinsResult.get(type_protein, 0) + 1
                door = "close"
                listPeptides[type_protein].append(pept)
    if type_protein in ["Globulin", "Triticin", "Serpin", "Avenin", "Hordein", "Secalin"]:
        if values in description.lower():
            if pept not in listPeptides[type_protein]:
                proteinsResult[type_protein] = proteinsResult.get(type_protein, 0) + 1
                door = "close"
                listPeptides[type_protein].append(pept)
    return(door)