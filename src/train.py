# -*- coding: utf-8 -*-
import csv
import pandas as pd
import re
import spacy
import sys

nlp = spacy.load('en_core_web_sm')

TRAINING_METADATA_FILE = "data/notes_metadata_trainingset_expanded.csv"
OUTPUT_FILE = "results/training_results.csv"

OLD_CHADSVASC_PATTERN = r"(?i)\bchads\b|\bchads2\b|\bchads-2\b|\bchads\b2\b|\bcha2\b"
NEW_CHADSVASC_PATTERN = r"(?i)\bcha2ds2vasc\b|\bcha2ds2-vasc\b|\bcha2ds2\svasc\b|\bchadsvasc\b|\bchads-vasc\b|\bchads2vasc\b|\bchads2-vasc\b|\bchads2\svasc\b|\bchads\svasc\b"


def train(metadata_frame, output_file):
	"""Extracts CHA₂DS₂-VASc scores from a patient note.

	Aggregates clinical notes by patient. For each patient, iterates through
	each phrase in each note, and tries to find target terms that indicate the
	term CHA₂DS₂-VASc, either in the earlier way of representation or in the
	current way of represenation. Extracts any numeral value within a span of
	15 characters or the end of phrase, whichever comes first, and labels it
	as the CHA₂DS₂-VASc score.

	Args:
		metadata_frame: The expanded metadata frame containing clinical notes
			for all mrns in the training set.
		output_file: The name of the file where the CHA₂DS₂-VASc scores will be
			written.

	Returns:
		None.

	Raises:
		None.
	"""

	patient_objs = []
	mrns = metadata_frame['mrn'].unique()
	for mrn in mrns:
		records = metadata_frame[metadata_frame['mrn'] == mrn]
		notes = []
		for row in records.itertuples():
			if not isinstance(row.text, str):
				continue
			notes.append((row.noteid, row.note_date, row.text))
		obj = {'mrn': mrn, 'notes': notes}
		patient_objs.append(obj)
	output_fieldnames = ['mrn', 'note_id', 'date', 'phrase', 'old targets', 'values', 'new targets', 'values']
	with open (output_file, 'w') as f:
		writer = csv.writer(f)
		writer.writerow(output_fieldnames)
		count = 0
		for patient in patient_objs:
			mrn = patient['mrn']
			for note in patient['notes']:
				note_id = note[0]
				date = note[1]
				note = note[2]
				doc = nlp(unicode(note, 'utf-8'))
				sentences = [str(i) for i in list(doc.sents)]
				for sentence in sentences:
					phrase = sentence.lower()
					old_targets = re.findall(OLD_CHADSVASC_PATTERN, sentence)
					old_target_positions = [m.span() for m in re.finditer(OLD_CHADSVASC_PATTERN, sentence)]
					old_chadsvasc_values = []
					for i in range(len(old_target_positions)-1):
						if old_target_positions[i][1]+15 > old_target_positions[i+1][0]:
							span = phrase[old_target_positions[i][1]:old_target_positions[i+1][0]]
						else:
							span = phrase[old_target_positions[i][1]:old_target_positions[i][1]+15]
						modifiers = re.findall(r"(?i)\bvasc\b", span)
						if len(modifiers) > 0:
							old_targets.remove(old_targets[i])
							continue
						old_chadsvasc_values = old_chadsvasc_values + re.findall(r"\d+", span)
					if len(old_target_positions) > 0:
						if len(phrase) > old_target_positions[-1][1]+15:
							last_span = phrase[old_target_positions[-1][1]:old_target_positions[-1][1]+15]
						else:
							last_span = phrase[old_target_positions[-1][1]:len(phrase)]
						modifiers = re.findall(r"(?i)vasc", last_span)
						if len(modifiers) > 0:
							old_targets.pop(-1)
						else:
							old_chadsvasc_values = old_chadsvasc_values + re.findall(r"\d+", last_span)
					new_targets = re.findall(NEW_CHADSVASC_PATTERN, sentence)
					new_target_positions = [m.span() for m in re.finditer(NEW_CHADSVASC_PATTERN, sentence)]
					new_chadsvasc_values = []
					for i in range(len(new_target_positions)-1):
						if new_target_positions[i][1]+15 > new_target_positions[i+1][0]:
							span = phrase[new_target_positions[i][1]:new_target_positions[i+1][0]]
						else:
							span = phrase[new_target_positions[i][1]:new_target_positions[i][1]+15]
						new_chadsvasc_values = new_chadsvasc_values + re.findall(r"\d+", span)
					if len(new_target_positions) > 0:
						if len(phrase) > new_target_positions[-1][1]+15:
							last_span = phrase[new_target_positions[-1][1]:new_target_positions[-1][1]+15]
						else:
							last_span = phrase[new_target_positions[-1][1]:len(phrase)]
						new_chadsvasc_values = new_chadsvasc_values + re.findall(r"\d+", last_span)
					if len(old_targets) > 0 or len(new_targets) > 0:
						if len(old_targets) > 0:
							old_target = old_targets[0]
						else:
							old_target = ''
						if len(old_chadsvasc_values) > 0:
							old_chadsvasc_value = old_chadsvasc_values[0]
						else:
							old_chadsvasc_value = ''
						if len(new_targets) > 0:
							new_target = new_targets[0]
						else:
							new_target = ''
						if len(new_chadsvasc_values) > 0:
							new_chadsvasc_value = new_chadsvasc_values[0]
						else:
							new_chadsvasc_value = ''
						line = [mrn, note_id, date, sentence, old_target, old_chadsvasc_value, new_target, new_chadsvasc_value]
						writer.writerow(line)
			count += 1
			sys.stdout.write('\rCompleted %d of %d, %d%%' % (count, len(patient_objs), count*100/len(patient_objs)))
			sys.stdout.flush()


def main():
	"""Main function.

	Reads training data. Extracts CHA₂DS₂-VASc scores. Writes the identifier
	terms and score value into a file.

	Usage:
	python train.py
	"""

	# read training data
	print 'Reading training metadata file...'
	metadata_frame = pd.read_csv(TRAINING_METADATA_FILE)

	# write CHA₂DS₂-VASc scores
	print 'Writing CHA₂DS₂-VASc scores...'
	train(metadata_frame=metadata_frame, 
			output_file=OUTPUT_FILE)


if __name__ == '__main__':
	main()