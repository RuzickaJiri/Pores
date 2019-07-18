import ast

with open("content.txt") as file:
    lines = file.readlines()
    str = lines[0]
print(str)

content = ast.literal_eval(str)
# print(content)
# print(list(content.values()))
# print(list(content.values())[0].values())
# print(list(list(content.values())[0].values())[0][2])
# print(list(content.keys()))

pores = []
for i in range(len(content)):
    if list(list(content.values())[i].values())[0][2] != 0:
        pores.append(list(content.keys())[i])

print(pores)
print(len(pores))
with open('pores.txt', 'w') as f:
    for item in pores:
        f.write("%s\n" % item)
