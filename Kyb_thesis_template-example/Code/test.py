# This is a comment
print("Hello World")

def num_input():
	return int(input("Enter number: "))

# Another comment
def validate_num(num):
	while num > 2:
		if num%2 == 0 and 2!=num:
			print("No prime")
			break
		else:
			print("Prime")
			break
	return num