# Testscript for Optimization Routines

# Set initial values, n is the dimension of the problem
n = 8
initial_values = 1.1 + 1:n / (4*n) 
#initial_values = 1:n

# Choose test function
test_function = 2
if (test_function == 1) { # see test_functions, squared squared Euclidean norm
  fname = sq_test
  grname = sq_grad_test
  hesname = sq_hesse_test
} else if (test_function == 2) { # see test_functions, Rosenbrock function
    fname = rosenbrock_advanced
    grname = grad_rosenbrock_advanced
} else if (test_function == 3) { # TODO, see proposal
  # still to implement
  print("TODO")
} else stop("No function specified")

if (1) {

  print("Start Gradient Descent")
  res <- grad_descent(start_val=initial_values, f_name=fname, grad_f_name = grname)

  print("Start Nonlinear CG")
  res2 <- conj_grad(start_val=initial_values, f_name=fname, grad_f_name = grname)

}

# TODO: Implement count for function and gradient evaluations (initialize counters as global variables)
