package lsq

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"os"
	"regexp"
)

// read data file utility

type DataFrame struct {
	Colnames []string
	Rows     [][]float64
	Ncol     int
	Nrow     int

	Rownames []string
}

func (df *DataFrame) String() string {
	s := fmt.Sprintf("DataFrame [%d x %d] = \n", df.Nrow, df.Ncol)
	digits := 3
	format := fmt.Sprintf(" %%.%df ", digits)

	var i, j int
	if len(df.Colnames) > 0 {
		for j = 0; j < df.Ncol; j++ {
			s += fmt.Sprintf(" %s ", df.Colnames[j])
		}
		s += "\n"
	}
	for i = 0; i < df.Nrow; i++ {
		for j = 0; j < df.Ncol; j++ {
			s += fmt.Sprintf(format, df.Rows[i][j])
		}
		s += "\n"
	}
	return s
}

/*
func main() {
	df, err := readData("fuelcons.dat")
	if err != nil {
		panic(err)
	}
	fmt.Printf("df = %v\n", df)
}
*/

func LineToFloatSlice(line string) []float64 {
	slc := LineToStringSlice(line)
	var z float64
	var a []float64
	for i := range slc {
		n, err := fmt.Sscanf(slc[i], "%f", &z)
		if n == 0 {
			a = append(a, math.NaN())
			continue
		}
		//fmt.Println("n", n)
		if err != nil {
			fmt.Printf("ERROR: %v\n", err)
			continue
		}
		a = append(a, z)
	}
	return a
}

var wordRegex = regexp.MustCompile(`(\S+)+`)

var firstWordRegex = regexp.MustCompile(`(\S+)`)

// right pad names with spaces, so all will be the same length
func NormalizeNameLengths(names []string) []string {
	maxLength := 0
	curlen := 0
	for i := range names {
		//fmt.Printf("names[%d] = '%v'   len(names[i]=%v   maxLength=%v\n", i, names[i], len(names[i]), maxLength)
		curlen = len(names[i])
		if curlen > maxLength {
			maxLength = curlen
		}
	}
	if maxLength > 0 {
		//fmt.Printf("maxLength = %v\n", maxLength)
		res := make([]string, len(names))
		format := fmt.Sprintf("%%-%ds", maxLength)
		for i := range names {
			res[i] = fmt.Sprintf(format, names[i])
			//fmt.Printf("res[i] = '%v'\n", res[i])
		}
		return res
	}
	return names
}

func LineToStringSlice(line string) []string {
	s := wordRegex.FindAllStringSubmatch(line, -1)

	//fmt.Printf("s = %v\n", s)

	slc := []string{}
	for i := range s {
		slc = append(slc, s[i][0])
	}
	//fmt.Printf("slc is %v\n", slc)
	return slc
}

func readData(fname string) (*DataFrame, error) {

	f, err := os.Open(fname)
	if err != nil {
		panic(fmt.Sprintf("error opening file '%s': %s", fname, err))
	}
	defer f.Close()
	df := &DataFrame{}
	r := bufio.NewReader(f)

	line, err := r.ReadString(10)
	if err != nil {
		panic(err)
	}
	df.Colnames = LineToStringSlice(line)
	df.Colnames = NormalizeNameLengths(df.Colnames)

	df.Ncol = int(len(df.Colnames))

	var rowname string
	for {
		line, err := r.ReadString(10) // 0x0A separator = newline
		if err == io.EOF {
			// do something here
			// fmt.Printf("line on eof: '%s'\n", line)
			break
		} else if err != nil {
			return nil, err
		}
		slice := LineToFloatSlice(line)
		//fmt.Printf("slice = %v\n", slice)
		df.Rows = append(df.Rows, slice)
		df.Nrow++

		// grab first word as rowname
		rn := wordRegex.FindStringSubmatch(line)
		if len(rn) > 0 {
			rowname = rn[0]
		} else {
			rowname = ""
		}
		df.Rownames = append(df.Rownames, rowname)
	}

	return df, nil
}
