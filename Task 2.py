import json
import re


def getting_test_sequences():
    sequences = []
    hla = []
    with open("test_texts.json", 'r') as json_file:
        data = json.load(json_file)
        for sequence in data:
            if sequence["hla"]:
                sequences.append(sequence["text"])
                hla.append(sequence["hla"])
    return sequences, hla


def data_preprocessing(text):
    text = re.sub(r'[^\w\s]', ' ', text)
    cleaned_text = text.replace("Human Leukocyte Antigen", "HLA" + " " * 20)

    for word in cleaned_text.split():
        if not word.isupper():
            cleaned_text = cleaned_text.replace(word, word[0].lower() + word[1:])

    for word in cleaned_text:
        if any(char.islower() for char in word or word.isalpha()):
            cleaned_text = cleaned_text.replace(word, '_' * len(word))

    return cleaned_text


def hla_finder_in_text(text):
    hla = []
    text = text.upper()
    matches = re.finditer('HLA', text)
    indexes = [match.start() for match in matches]

    for index in indexes:
        i = 3
        while not text[index + i].isalpha() and index + i <= len(text) - 4:
            i += 1

        if index + i not in indexes:
            hla.append({"gene": "", "allele": None, "protein": None, "positions": [[index + changes]]})

            while text[index + i].isalpha() or (text[index + i].isdigit() and len(hla[-1]["gene"]) > 1):
                hla[-1]["gene"] += text[index + i]
                i += 1

            if hla[-1]["gene"] == '' or hla[-1]["gene"][0] not in ('A', 'B', 'C', 'D', 'E'):
                hla.pop()
                continue

            positions_end = index + i
            skip = False
            tries = 0
            while not text[index + i].isdigit() and index + i <= len(text) - 2:
                if index + i in indexes:
                    skip = True
                if text[index + i] in 'DQB':
                    tries += 1
                i += 1
            if skip or tries >= 2:
                hla[-1]["positions"][-1].append(positions_end)
                continue

            allele = ''
            while text[index + i].isdigit() and len(allele) < 2 and text[index + i]:
                allele += text[index + i]
                i += 1

            if allele != '':
                hla[-1]["allele"] = allele
                positions_end = index + i

            while not text[index + i].isdigit() and index + i <= len(text) - 2:
                i += 1
            protein = ''

            while text[index + i].isdigit() and len(protein) < 3:
                protein += text[index + i]
                i += 1

            if protein != '':
                hla[-1]["protein"] = protein
                hla[-1]["positions"][-1].append(index + i)
            else:
                hla[-1]["positions"][-1].append(positions_end)
        else:
            continue

    return hla


def many_genes_without_hla(text):
    global changes
    changes = 0
    text_for_check = text.replace("_", " ")
    text_for_check = text_for_check.split()
    is_hla = False
    for i in range(len(text_for_check)):
        if text_for_check[i] == "HLA":
            j = text.find("HLA")
            is_hla = True
            matches = re.finditer(r'\b{}\b'.format(re.escape(text_for_check[i+1])), text)
            index = [match.start() for match in matches]
            text_list = list(text)
            text_list[index[0]:index[0] + len(text_for_check[i+1])] = " " * len(text_for_check[i+1])
            text = ''.join(text_list)
            break

    if is_hla and i + 2 < len(text_for_check):
        new_hla_index = text.find(text_for_check[i+2])
        text = (text[:j] + text[new_hla_index-3:new_hla_index] + text[j+3:new_hla_index-4] + text[j:j+4] +
                text[new_hla_index:])
        changes = 4

    for i in range(len(text_for_check)):
        if text_for_check[i].isdigit():
            index = text.find(text_for_check[i])
            text_list = list(text)
            text_list[index:index + len(text_for_check[i])] = " " * len(text_for_check[i])
            text = ''.join(text_list)
        else:
            break

    return text


sequences, hla_test = getting_test_sequences()
count = 0

for index, sequence in enumerate(sequences):
    global changes
    changes = 0
    cleaned_sequence = data_preprocessing(sequence)
    hla_train = hla_finder_in_text(cleaned_sequence)
    if hla_train == hla_test[index]:
        count += 1
    else:
        new_cleaned_sequence = many_genes_without_hla(cleaned_sequence)
        new_hla_train = hla_finder_in_text(new_cleaned_sequence)
        hla_train.extend(new_hla_train)
        while hla_train != hla_test[index]:
            new_cleaned_sequence = many_genes_without_hla(new_cleaned_sequence)
            if new_cleaned_sequence == sequence:
                print(new_cleaned_sequence)
                print(hla_train, hla_test[index], sep='\n', end='\nNEW\n')
                break
            new_hla_train = hla_finder_in_text(new_cleaned_sequence)
            hla_train.extend(new_hla_train)
        else:
            count += 1

print(count,"right identified sequences out of", len(sequences))
