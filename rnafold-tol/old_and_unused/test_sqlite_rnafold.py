import sqlite_rnafold as db

db_, connection = db.get_db("test")

s1 = db.Store("test")
print(s1.get_float("X12345.2", 101))

try:
    s1.set_float("X12345.2", 101, 1.01010101)
except db.IntegrityError as e:
    print(e)

print(s1.get_float("X12345.2", 101))




