function rowSwap(A,row1,row2)
	A[row1,:], A[row2,:] = A[row2,:], A[row1,:]
end
