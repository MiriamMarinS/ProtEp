#For epitopes with original sequence. Method is the type of epitope: with original sequence, deaminated1 or deaminated2 sequence, p31 sequence, wage sequence, moAb sequence or ias sequence.
from Bio import SeqIO

class Epitope():
    def __init__(self, epitope_file, method):
        self.epitope_file = epitope_file
        self.method = method
        self.library = self.epitope_library()

    def epitope_library(self):
        self.epitope_sequences = {}
        type_epitopes = [("a", "alpha"), ("g", "gamma"), ("w", "omega"), ("L", "LMW"), ("H", "HMW"), ("hor", "hordein"), ("sec", "secalin"), ("ave", "avenin"), ("p31", "p31"), ("a1", "a1"), ("g12", "g12"), ("r5", "r5")]
        for record in SeqIO.parse(self.epitope_file, "fasta"):
            if str(record.id).split(";")[1] == self.method:
                if self.method == "original" or self.method == "deaminated1" or self.method == "deaminated2":
                    if "glia" in str(record.id) or "glut" in str(record.id):
                        epitope = [type_epitope[1] for type_epitope in type_epitopes if str(record.id).split(";")[0].split("_")[2].startswith(type_epitope[0])][0]
                    elif any(i in str(record.id) for i in ["hor", "sec", "ave"]):
                        epitope = [type_epitope[1] for type_epitope in type_epitopes if str(record.id).split(";")[0].split("_")[1].startswith(type_epitope[0])][0]
                    else:
                        epitope = [type_epitope[1] for type_epitope in type_epitopes if str(record.id).split(";")[0].split("_")[0].startswith(type_epitope[0])][0]

                elif self.method == "moAb":
                    epitope = [type_epitope[1] for type_epitope in type_epitopes if str(record.id).split(";")[0].split("_")[0].strip() == type_epitope[0]][0]
                else:
                    epitope = "non_classified"
                self.epitope_sequences[str(record.id).split(";")[0]] = (str(record.seq), epitope)

        return(self.epitope_sequences)

    def epitope_library_remove_duplication(self):
        # In the case of original, deaminated1, deaminated2 and moAb epitopes methods, the epitopes are grouped in classes.
        # If epitopes in the same class have the same sequences, they will be considered as one epitope. In other methods, the duplication will be removed in the group of all epitopes, because they are not grouped in classes.

        #Find duplicated sequences in epitope library.
        duplicated_sequences = {}
        for epitope_id, value in self.library.items():
            epitope_seq = value[0].strip()
            epitope_type = value[1]
            duplicated_sequences.setdefault((epitope_seq, epitope_type), []).append(epitope_id)

        #New dictionary with epitopes grouped by their sequences separated in original, deaminated1 and deaminated2 groups, or unclassified group for the other methods.
        self.new_epitope_library = {}
        for k, v in duplicated_sequences.items():
            self.new_epitope_library[';'.join(v)] = k


        #The name of epitopes with duplicated sequences are concatenated by ;. But the name is too long. So we kept the first epitope name, and it can be checked in epitope library printed.
        self.epitope_sequences2 = {}
        for k, v in self.new_epitope_library.items():
            if len(k.split(";")) > 1:
                self.epitope_sequences2[k.split(";")[0]] = v
            else:
                self.epitope_sequences2[k] = v


        return(self.epitope_sequences2)