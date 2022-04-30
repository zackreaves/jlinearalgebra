# Function that uses gaussian elimination to simplify matrices, all inputs must be in the form of a numerical matrix(A).
function gaussElim(A)

	# Parameters
	rows, columns = size(A)
	if typeof(A) == Matrix{Int}
		A = A//1
	else
		A = convert(Matrix{Float64}, A)
	end

	# Loop that should ideally swap rows as needed, then reduce pivots to one, then eliminate all other non-zero numbers in the column.

	for column = 1:columns
		for row = column:rows

			# Sorting Loop

			swapRow = row
			while A[row,column] == 0 && swapRow < rows
				swapRow += 1
				A[row,:], A[swapRow,:] = A[swapRow,:], A[row,:]
				if swapRow == rows && A[row,column] == 0 && column != columns
					swapRow = row
					column += 1
				end
			end

			# Reducing pivots to one

			if A[row,column] != 0
				A[row,:] /= A[row,column]
			#elseif sum(A[row,:]) != 0
			#	while A[row,column] == 0
			#		column += 1
			#	end
			#	A[row,:] /= A[row,column]
			end

			# Zeroing out columns

			for zeRow = 1:rows
				if zeRow != row && sum(A[row,:]) != 0
					A[zeRow,:] -= A[row,:] * A[zeRow,column]
				end
			end
		end
	end

	if typeof(A) == Matrix{Rational{Int}} && sum(mod.(numerator.(A), denominator.(A))) == 0
		A = convert(Matrix{Int}, A)
	else
		A = convert(Matrix{Any}, A)
		for row = 1:rows
			for column = 1:columns
				if typeof(A[row,column]) == Rational{Int} && mod(numerator(A[row,column]), denominator(A[row,column])) == 0
					A[row,column] = convert(Int, A[row,column])
				end
			end
		end
	end

	return A
end
