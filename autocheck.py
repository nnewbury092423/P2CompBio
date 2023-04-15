from testing.tests_autocheck import *
from os.path import exists

def check_identity():
    student_info = {'Name':None,'PID':None,'Github ID':None}
    with open("student_info.txt",'r') as fin:
        for line in fin:
            x,y = line.strip().split(':')
            if x in student_info and y:
                student_info[x.strip()] = y.strip()
    for key in student_info:
        assert student_info[key] is not None, "Invalid submission: missing " + key + " in student_info.txt"
    assert exists("writeup_" + student_info["PID"] + ".pdf"), "Missing writeup file. Make sure you name the file as writeup_<yourPID>.pdf"
    return student_info

if __name__ == "__main__":
    #student_info = check_identity()
    #print("Running autocheck for student " + student_info["PID"])
    unittest.main()
