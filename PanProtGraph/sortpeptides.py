import csv
import sys

def load_file(filename):
  """Loads a file with two columns of strings separated by ';'."""
  with open(filename, 'r') as f:
    reader = csv.reader(f, delimiter=';')
    data = []
    for row in reader:
      data.append(row)
  return data

def sort_rows(data):
  """Sorts the rows based on the 1st and 2nd strings."""
  sorted_data = []
  for row in data:
    sorted_data.append((row[0], row[1]))
  sorted_data.sort()
  return sorted_data

def main():
  """Loads the file, sorts the rows, and removes duplicated pairs of strings."""
  data = load_file(sys.argv[1])
  print("File loaded")
  sorted_data = sort_rows(data)
  print("sorted")
  print(len(sorted_data))

  unique_data = []
  prevrow = [[], []]
  for row in sorted_data:
      if row[0] != prevrow[0] or row[1] != prevrow[1]:
        unique_data.append(row)
        prevrow = row
  print(len(unique_data))

  with open("uniquepeptides_by_ethnics.txt", "w") as f:
    for row in unique_data:
        f.write(f"{row[0]};{row[1]}\n")

if __name__ == '__main__':
  main()
