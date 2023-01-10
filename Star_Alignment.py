import numpy as np

sequences = []
array_length = int(input())

for i in range(0, array_length):
  temp = input()
  sequences.append(temp)

def global_align(x, y, s_match = 3, s_mismatch = -1, s_gap = -2):
    A = []
    for i in range(len(y) + 1):
        A.append([0] * (len(x) + 1))
    for i in range(len(y) + 1):
        A[i][0] = s_gap * i
    for i in range(len(x) + 1):
        A[0][i] = s_gap * i
    for i in range(1, len(y) + 1):
        for j in range(1, len(x) + 1):
            A[i][j] = max(
                A[i][j - 1] + s_gap,
                A[i - 1][j] + s_gap,
                A[i - 1][j - 1] + (s_match if (y[i - 1] == x[j - 1] and y[i - 1] != '-') else 0) + (
                    s_mismatch if (y[i - 1] != x[j - 1] and y[i - 1] != '-' and x[j - 1] != '-') else 0) + (
                    s_gap if (y[i - 1] == '-' or x[j - 1] == '-') else 0)
            )
    align_X = ""
    align_Y = ""
    i = len(x)
    j = len(y)
    while i > 0 or j > 0:
        current_score = A[j][i]
        if i > 0 and j > 0 and (
                ((x[i - 1] == y[j - 1] and y[j - 1] != '-') and current_score == A[j - 1][i - 1] + s_match) or
                ((y[j - 1] != x[i - 1] and y[j - 1] != '-' and x[i - 1] != '-') and current_score == A[j - 1][
                    i - 1] + s_mismatch) or
                ((y[j - 1] == '-' or x[i - 1] == '-') and current_score == A[j - 1][i - 1] + s_gap)
        ):
            align_X = x[i - 1] + align_X
            align_Y = y[j - 1] + align_Y
            i = i - 1
            j = j - 1
        elif i > 0 and (current_score == A[j][i - 1] + s_gap):
            align_X = x[i - 1] + align_X
            align_Y = "-" + align_Y
            i = i - 1
        else:
            align_X = "-" + align_X
            align_Y = y[j - 1] + align_Y
            j = j - 1
    return (align_X, align_Y, A[len(y)][len(x)] )

def score(sequences):
  score_value = int(0)
  for i in range(0, len(sequences)-1):
    for j in range(i+1, len(sequences)):
      for k in range(0, len(sequences[0])):
        if (sequences[i][k] == "-") and (sequences[j][k] == "-"):
          score_value += 0
        elif (sequences[i][k] == "-") or (sequences[j][k] == "-"):
          score_value += -2
        elif sequences[i][k] == sequences[j][k] :
          score_value += 3
        else :
          score_value += -1
  return int(score_value)

def dash_placed(s1, s2, s3, k=0):
  pos = []
  i= 0
  while i < len(s1):
    if s2[i] != s1[i] :
      pos.append(i)
      s2 = s2[0:i] + s2[i+1:]
    else :
      i += 1

  dash_at_last = len(s2) - len(s1)

  pos = np.flip(pos)

  for i in range(0, len(pos)):
    s3 = s3[0:pos[i]] + '-' + s3[pos[i]:]

  for _ in range(0,dash_at_last):
    if k == 1 :
      while True:
        pass
    s3 = s3 + "-"

  return s3

def MSA(sequences):

  similarity_matrix = np.zeros( (len(sequences), len(sequences)) )

  for i in range(0, len(sequences)):
    for j in range(0, len(sequences)):
      if i != j :
        similarity_matrix[i, j] = global_align(sequences[i], sequences[j])[2]
  sum_of_similarities = np.zeros((1, len(sequences)))

  for i in range(0, len(sequences)):
    sum_of_similarities[0][i] = np.sum(similarity_matrix[i])

  high_to_low = np.flip(sum_of_similarities.argsort())[0]
  center_sequence_number = high_to_low[0]
 
  high_to_low = []
  center_sequence = sequences[center_sequence_number]

  temp = similarity_matrix[center_sequence_number,:].argsort()
  temp = np.flip(temp)
  for i in range(len(sequences)):
    if temp[i] != center_sequence_number:
      high_to_low.append(temp[i])

  for _ in range(10):
    for i in range(len(high_to_low)-1):
      if (similarity_matrix[center_sequence_number, high_to_low[i]] == similarity_matrix[center_sequence_number, high_to_low[i+1]]):
        if(high_to_low[i] > high_to_low[i+1]):
          temp = high_to_low[i]
          high_to_low[i] = high_to_low[i+1]
          high_to_low[i+1] = temp

  final = []
  for i in range(0, array_length):
    final.append("")


  final[center_sequence_number] = global_align(sequences[high_to_low[0]], center_sequence)[1]
  final[high_to_low[0]] = global_align(sequences[high_to_low[0]], center_sequence)[0]

  for i in high_to_low[1:] :
    temp = final[center_sequence_number]
    final[center_sequence_number] = global_align(sequences[i], final[center_sequence_number])[1]
    
    for j in range(0, len(sequences)):
      if (final[j] != "") and (j != center_sequence_number) :
        final[j] = dash_placed(temp, final[center_sequence_number], final[j])

    final[i] = global_align(sequences[i], final[center_sequence_number])[0]
  
  return final

def Block_MSA(final):
  new_sequences = []
  for _ in range(0, array_length):
      new_sequences.append("")

  temp_final = []
  for _ in range(0, array_length):
      temp_final.append("")

  same_columns = []

  i= 0
  while i < len(final[0]):
    if_all_seq_are_same = 1
    for j in range(1, len(sequences)):
      if (final[0][i] == final[j][i]) and (final[0][i] != "-"):
        if_all_seq_are_same += 1
        if if_all_seq_are_same == len(sequences):
          same_columns.append(i)
        if (len(same_columns) == 2) and (same_columns[1] > (same_columns[0]+2)) :
          for m in range(0, len(sequences)) :
            for n in range(same_columns[0]+1, same_columns[1]) :
              if final[m][n] != "-" :
                new_sequences[m] = new_sequences[m] + final[m][n]

          for p in range(0, len(sequences)):
            temp_final[p] = final[p][0:same_columns[0]+1] + MSA(new_sequences)[p] + final[p][same_columns[1]:]
          if score(temp_final) > score(final) :
            final = temp_final
            i = -1
          else :
            i -= 1         
          temp_final = []
          for _ in range(0, array_length):
            temp_final.append("")
          new_sequences = []
          for _ in range(0, array_length):
            new_sequences.append("")
          same_columns = []
        elif (len(same_columns) == 2) and (same_columns[1] <= (same_columns[0]+2)) :
          same_columns = []
    i += 1

  new_sequences = []
  for _ in range(0, array_length):
      new_sequences.append("")

  temp_final = []
  for _ in range(0, array_length):
      temp_final.append("")

  same_columns = []

  i= 0
  while i < len(final[0]):
    if_all_seq_are_same = 1
    for j in range(1, len(sequences)):
      if (final[0][i] == final[j][i]) and (final[0][i] != "-"):
        if_all_seq_are_same += 1
        if if_all_seq_are_same == len(sequences):
          same_columns.append(i)
    i += 1

  if (len(same_columns) >= 2) :
    if (same_columns[0] >= 2) :
      for m in range(0, len(sequences)) :
        for n in range(0, same_columns[0]) :
          if final[m][n] != "-" :
            new_sequences[m] = new_sequences[m] + final[m][n]
      for p in range(0, len(sequences)):
        temp_final[p] = MSA(new_sequences)[p] + final[p][same_columns[0]:]
      if score(temp_final) > score(final) :
        final = temp_final

    new_sequences = []
    for _ in range(0, array_length):
      new_sequences.append("")

    temp_final = []
    for _ in range(0, array_length):
      temp_final.append("")

    if (same_columns[-1] < len(final[0]) - 2) :
      for m in range(0, len(sequences)) :
        for n in range(same_columns[-1]+1, len(final[0]) ) :
          if final[m][n] != "-" :
            new_sequences[m] = new_sequences[m] + final[m][n]
      for p in range(0, len(sequences)):
        temp_final[p] = final[p][:same_columns[-1]+1] + MSA(new_sequences)[p]
      if score(temp_final) > score(final) :
        final = temp_final

  elif (len(same_columns) == 1) : 
    if (same_columns[0] >= 2) :
      for m in range(0, len(sequences)) :
        for n in range(0, same_columns[0]) :
          if final[m][n] != "-" :
            new_sequences[m] = new_sequences[m] + final[m][n]
      for p in range(0, len(sequences)):
        temp_final[p] = MSA(new_sequences)[p] + final[p][same_columns[0]:]
      if score(temp_final) > score(final) :
        final = temp_final

    new_sequences = []
    for _ in range(0, array_length):
      new_sequences.append("")

    temp_final = []
    for _ in range(0, array_length):
      temp_final.append("")

    if (same_columns[0] < len(final[0]) - 2) :
      for m in range(0, len(sequences)) :
        for n in range(same_columns[0]+1, len(final[0]) ) :
          if final[m][n] != "-" :
            new_sequences[m] = new_sequences[m] + final[m][n]
      for p in range(0, len(sequences)):
        temp_final[p] = final[p][:same_columns[0]+1] + MSA(new_sequences)[p]
      if score(temp_final) > score(final) :
        final = temp_final

  return final

final = MSA(sequences)

while True :
  temp = final
  final = Block_MSA(final)
  if temp == final :
    break

print(score(final))

for i in range(len(sequences)):
  print(final[i])