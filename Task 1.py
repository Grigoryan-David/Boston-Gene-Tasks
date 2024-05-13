import yaml
import re
import json


def genes_from_yaml():
    with open("genes.yaml", 'r') as yaml_file:
        data = yaml.safe_load(yaml_file)

    for i in data:
        for j, synonym in enumerate(i['synonyms']):
            i['synonyms'][j] = data_preprocessing(synonym)
        i['synonyms'] = list(set(i['synonyms']))

    return data


def genes_finder_in_text(text):
    genes_in_sequence = {"genes": []}
    text = text.replace("can", "   ")
    text = text.lower()
    for gene in genes_data:
        already_is = False
        for synonym in gene["synonyms"]:
            length = len(synonym)
            synonym = r'\b{}\b'.format(synonym)
            matches = re.finditer(synonym.lower(), text)
            indexes = [match.start() for match in matches]
            if indexes:
                if not already_is:
                    genes_in_sequence["genes"].append({"name": gene["name"], "positions": []})
                    already_is = True
                    for index in indexes:
                        genes_in_sequence["genes"][-1]["positions"].append([index, index + length])
                else:
                    for index in indexes:
                        genes_in_sequence["genes"][-1]["positions"].append([index, index + length])

    for gene in genes_in_sequence["genes"]:
        gene["positions"].sort()
        for position1 in gene["positions"]:
            for position2 in gene["positions"]:
                if position1[0] == position2[0] and position1[1] > position2[1]:
                    gene["positions"].remove(position2)

    return genes_in_sequence


def getting_test_sequences():
    sequences = []
    genes = []
    with open("test_texts.json", 'r') as json_file:
        data = json.load(json_file)
        for sequence in data:
            if sequence["genes"]:
                cleaned_sequence = data_preprocessing(sequence["text"])
                sequences.append(cleaned_sequence)
                genes.append(sequence["genes"])
    return sequences, genes


def data_preprocessing(text):
    pattern = r'[^\w\s]'
    cleaned_text = re.sub(pattern, ' ', text)
    return cleaned_text


genes_data = genes_from_yaml()
test_sequences, test_genes = getting_test_sequences()
count = 0
for index, sequence in enumerate(test_sequences):
    genes_in_text = genes_finder_in_text(sequence)
    if genes_in_text["genes"] == test_genes[index]:
        count += 1
    else:
        print(sequence, genes_in_text["genes"], test_genes[index], sep="\n")

print(count,"right identified sequences out of", len(test_sequences))
