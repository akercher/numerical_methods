def main():
	import sys
	filename = sys.argv[1]
	f = open(filename, 'r')
	print 'reading in ', filename
	without_whitespace = '\n'.join([line.strip() for line in f.readlines()])
	
	f = open(filename, 'w')
	f.write(without_whitespace)
	print 'writing back to ', filename, ' without leading and trailing whitespace'

if __name__ == '__main__':
	main()
