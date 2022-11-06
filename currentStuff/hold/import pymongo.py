import pymongo
import mysql.connector

print ("version:", pymongo.version)



mydb = mysql.connector.connect(
  host="localhost",
  user="rna",
  password="rna"
)

print(mydb)