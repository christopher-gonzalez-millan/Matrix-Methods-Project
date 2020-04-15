from polynomial import *
import fieldmath


class GeneralizedReedSolomon(object):
    """
    This class is an implementation of the Generalized Reed Solomon error correcting code as presented by Yehuda Lindell
    in section 6 of the following lecture notes: http://u.cs.biu.ac.il/~lindell/89-662/main-89-662.html
    Reed-Solomon is a linear [n,k,n-k+1]-code meaning we can correct up to (n-k)//2 errors.
    """
    generator_matrix = None
    parity_check_matrix = None
    vp_arr = None  # Note: this is used for generation while v_arr is for parity checks

    def __init__(self, f, k, n=None, conventional_creation=False, alpha=None, v_arr=1, print_matrices=False):
        """
        With the provided parameters we can generate the generator and parity check matrices and generator polynomials.
        It is up to the user to ensure all provided inputs for __init__ encode and decode satisfy the field.
        Exceptions will be raised only by field class.
        :param f: Provided field class to perform all mathematical operations under. The recommended field is
        BinaryField(0x11d). Any valid field is allowed. Note neither this class nor Field check to ensure a valid field
        :param alpha: The code locators. The size of this array determines what n is. These values are entirely
        arbitrary, but must exist in the field f. a single value for alpha should be used if conventional_creation=True
        conventional_creation is explaining in section 6.2 under Classical Reed-Solomon.
        :param n: This parameter will automatically generate the alpha array, to contain n elements. It is recommended
        that conventional_creation=True when using this option. Note: either n or alpha must be not None. If both are
        provided alpha will be selected as default.
        :param conventional_creation: if True both a single alpha and n must be provided.
        :param v_arr: Column multipliers for H. Can either be a list of non-zero elements or a single element of the
        field. Ex v_arr=1 to use a normalized GRS, where 1 is the 1 element of the provided field.
        :param print_matrices: boolean on whether code should print G and H once they are created
        """
        # First we need to make sure our field is valid.
        if isinstance(f, fieldmath.Field):
            self.f = f
        else:
            raise Exception("Please provide a valid Field")

        # Next we want to work on our alpha array
        self.alpha_arr = []
        if alpha and isinstance(alpha, list):
            self.alpha_arr = alpha  # This assumes the user knows that they're doing
        elif n is not None:
            if conventional_creation:
                if alpha and not isinstance(alpha, list):
                    for i in range(1, n+1):
                        self.alpha_arr.append(fieldmath.pow_over_field(alpha, i, self.f))
                else:
                    raise Exception("alpha seems to be incorrect, ensure it is a single element of the specified field")
            else:
                try:
                    self.alpha_arr = [self.f.multiply(self.f.one(), i) for i in range(2, n + 2)]
                except:
                    raise Exception("Failed to create alpha array. Try conventional_creation=True")
        else:
            raise Exception("Either the parameter alpha_arr or n must be included")
        print(self.alpha_arr)
        # Set basic parameters for the code based on other inputs
        self.n = len(self.alpha_arr)
        self.k = k
        self.d = self.n - self.k + 1

        # save column multipliers for the parity check matrix
        if isinstance(v_arr, int):  # TODO: change from type int to anything really
            self.v_arr = [self.f.multiply(v_arr, self.f.one())]*self.n
        elif isinstance(v_arr, list):
            self.v_arr = v_arr  # Assumes user knows what they're doing
        else:
            raise TypeError("v_arr must be of type int or list")

        # A few checks to make sure everything is correct
        # q >= n not checked. Assumes user knows what they're doing
        if not self.k <= self.n:
            raise Exception(f"Does not satisfy k <= n: k={self.k}, n={self.n}")

        # Now we can create the matrices
        print(f"Initializing linear GRS [{self.n}, {self.k}, {self.d}]-code.")
        print(f"This GRS code can correct up to {int((self.d-1)/2)} errors.\n")
        self.create_matrices()
        if print_matrices:
            print(f"Generator matrix: \n{self.generator_matrix}")
            print(f"Parity Check matrix: \n{self.parity_check_matrix}\n")

    def encode(self, msg, use_poly=True):
        """
        There are two methods of doing this: use generator matrix or generator polynomial.
        Default is the polynomial, both work and result in the same outcome. Although matrix multiplication is rather
        simple for often than not RS is defined with the generator polynomial method.
        :param msg: input message of size k. Elements must be elements of the field f.
        :param use_poly: boolean to specify the use of the generator polynomial over generator matrix
        :return: encoded message of size n
        """
        if len(msg) != self.k:
            raise Exception(f"Input message must be of size k, k={self.k}")

        if use_poly:
            msg_poly = Polynomial(*msg, f=self.f)
            encoded_msg = []
            for i in range(0, self.n):
                encoded_msg.append(self.f.multiply(self.vp_arr[i], msg_poly(self.alpha_arr[i])))
            return encoded_msg
        else:
            # m*G where m is our msg
            msg_matrix = fieldmath.Matrix(1, self.k, self.f)
            for i in range(self.k):
                msg_matrix.set(0, i, msg[i])

            return msg_matrix * self.generator_matrix

    def decode(self, msg):
        """
        Decodes and error corrects. This is by far the most difficult part. This implementation uses the
        Peterson-Gorenstein-Zierier GRS Decoding algorithm as presented in section 6.3 and 6.3.4 of  the lecture notes.
        All calculations should be done in provided field. Very generally this approach solves the key equations (6.3.3)
        by find the kernel of the matrix created from the last d - 1 - tau equations. these will represent the
        coefficients of the lambda polynomial. We can then plug these back into the other key equations and fine the
        coefficients for the gamma polynomial. Note that these polynomials do not represent the actualy polynomials we
        need and to find the polynomials we need we can divide lambda by the gcd of lambda and gamma. This will van then
        be used to find where errors are and we can then solve our syndrome equations to find what the errors are and
        we then correct them and have our corrected encoded message.
        :param msg: Input message of size n
        :return: returns decoded and error corrected message of size k
        """
        # First find syndrome (S_0, ... , S_d-2)
        msg_syndrome = self.syndrome(msg)
        msg_matrix = fieldmath.create_matrix([msg], self.f)

        if msg_syndrome.any():
            tau = (self.d - 1) // 2
            syndrome_matrix = fieldmath.Matrix(self.d - 1, tau + 1, self.f, zeros=True)
            for i in range(self.d - 1):
                for j in range(i, max(-1, i - tau - 1), -1):
                    syndrome_matrix.set(i, i - j, msg_syndrome.get(0, j))

            # With the syndrome matrix we want to solve the last d-1-tau equations
            # To find lambda_poly
            lam_eqs = syndrome_matrix.get_sub_matrix(tau, None, None, None)
            lam_kernel_space = lam_eqs.kernel_space()
            if lam_kernel_space != 0:
                lam_coeff_matrix = lam_kernel_space * fieldmath.create_matrix([[1]] * lam_kernel_space.column_count(),
                                                                              self.f)
                lam_coeff = lam_coeff_matrix.to_list(single=True)
                lambda_poly = Polynomial(*lam_coeff, f=self.f)

                # plug lambda back into key equations to find gamma_poly
                gamma_coeff_matrix = syndrome_matrix.get_sub_matrix(None, tau, None, None)*lam_coeff_matrix
                gamma_poly = Polynomial(*gamma_coeff_matrix.to_list(single=True), f=self.f)

                # Calculate the GCD
                gcd_lg = poly_gcd(lambda_poly, gamma_poly)

                # divide to find error_locator_poly
                error_locator_poly, rem = divmod(lambda_poly, gcd_lg)

                # TODO: account for this
                if rem != Polynomial(self.f.zero(), f=self.f):
                    raise Exception("Could not decode message. Not GCD found.")

                # Find where the errors are by finding when Lambda(alpha_arr_k^-1) == 0
                error_locations = []
                for j in range(len(self.alpha_arr)):
                    alpha_inv = self.f.reciprocal(self.alpha_arr[j])
                    if error_locator_poly(alpha_inv) == 0:
                        error_locations.append(j)

                if len(error_locations) != 0:
                    # Now to actually figure out what the errors are. Solve Syndrome for e_j for each j in error
                    # locations
                    err_matrix = fieldmath.Matrix(self.d - 1, len(error_locations), self.f)
                    for r in range(err_matrix.row_count()):
                        for c in range(err_matrix.column_count()):
                            val = self.f.multiply(self.v_arr[error_locations[c]],
                                                  fieldmath.pow_over_field(self.alpha_arr[error_locations[c]], r,
                                                                           self.f))
                            err_matrix.set(r, c, val)

                    # This next line has resulted in an error once in multiple thousands of test and I haven't figured
                    # why. And I have not be able to reproduce it.
                    try:
                        errors = fieldmath.solve_ax_b(err_matrix, msg_syndrome.transpose())
                    except Exception as e:
                        print(f"Could not solve and find errors using solve_ax_b: {e}")
                        errors = fieldmath.create_matrix([[0]]*err_matrix.column_count(), self.f)

                    # finally fix the errors
                    for i in range(len(error_locations)):
                        msg_matrix.set(0, error_locations[i],
                                       self.f.subtract(msg_matrix.get(0, error_locations[i]), errors.get(i, 0)))

        # Last step: give provided msg that created the encoded msg. Because we aren't using G and H in standard form
        # We need to solve for what msg provided us with the output.
        # Ive noticed some weird behavior when the number of errors is greater than the number of fixable errors.
        if not self.syndrome(msg_matrix).any():
            return fieldmath.solve_ax_b(self.generator_matrix.transpose(), msg_matrix.transpose()).to_list(single=True)
        else:
            try:
                return fieldmath.solve_lstsq(self.generator_matrix.transpose(), msg_matrix.transpose()).to_list(single=True)
            except:
                return [0]*self.k

    def syndrome(self, msg):
        """
        Computes the syndrome of a code. This is used in the decode function.
        :param msg:
        :return: S_H(msg) = (S_0, ... , S_{d-2})
        """
        if isinstance(msg, list):
            if len(msg) != self.n:
                raise Exception("Input message must be of size n.")
            msg_matrix = fieldmath.Matrix(1, self.n, self.f)
            for i in range(self.n):
                msg_matrix.set(0, i, msg[i])
        elif isinstance(msg, fieldmath.Matrix):
            if msg.column_count() != self.n:
                raise Exception("Input message must be of size n.")
            msg_matrix = msg.copy()
        else:
            raise TypeError()

        # y*H^T where y is our msg
        syndrome = msg_matrix * self.parity_check_matrix.transpose()
        return syndrome

    def create_matrices(self):
        """
        This acts as a helper function to create the generator and parity check matrices
        TODO: Added option for matrices to be in standard form G = (I|X), H = (X^T|I)
        """
        self.generator_matrix = fieldmath.Matrix(self.k, self.n, self.f)
        self.parity_check_matrix = fieldmath.Matrix(self.n - self.k, self.n, self.f)

        # Create Parity Check matrix
        self.parity_check_matrix = self.create_parity_check_matrix(self.k)

        k_one_t = self.create_parity_check_matrix(1)
        k_one_kernel_space = k_one_t.kernel_space()
        if isinstance(k_one_kernel_space, int):
            raise Exception("Kernel space = 0. Could not find generator matrix. Try different alpha.")

        self.vp_arr = (k_one_kernel_space * fieldmath.create_matrix([[1]] * k_one_kernel_space.column_count(),
                                                                    self.f)).to_list(single=True)

        # Create Generator matrix from alpha and vp_arr
        for i in range(self.generator_matrix.row_count()):  # rows
            for j in range(self.generator_matrix.column_count()):  # columns
                val = self.f.multiply(self.vp_arr[j], fieldmath.pow_over_field(self.alpha_arr[j], i, self.f))
                self.generator_matrix.set(i, j, val)

        # Acts as a final check that the generator and parity check matrices are valid
        if (self.parity_check_matrix*self.generator_matrix.transpose()).any():
            raise Exception("Generator and Parity Check matrices do no align.")
        if not isinstance(self.generator_matrix.transpose().kernel_space(), int):
            raise Exception("Generator matrix has LD rows.")
        if not isinstance(self.parity_check_matrix.transpose().kernel_space(), int):
            raise Exception("Generator matrix has LD rows.")

    def create_parity_check_matrix(self, k):
        """
        Helper function to create parity check matrix for a give k
        :param k:
        :return: parity check matrix
        """
        pc_matrix = fieldmath.Matrix(self.n - k, self.n, self.f)

        for i in range(0, pc_matrix.row_count()):  # rows
            for j in range(0, pc_matrix.column_count()):  # columns
                val = self.f.multiply(self.v_arr[j], fieldmath.pow_over_field(self.alpha_arr[j], i, self.f))
                pc_matrix.set(i, j, val)
        return pc_matrix


if __name__ == '__main__':
    # The below is an example implementation of a GRS [15, 10, 6]-code in the binary field of 256 elements.
    field = fieldmath.BinaryField(0x11d)
    grs = GeneralizedReedSolomon(f=field, k=10, n=15, v_arr=1, alpha=0x2, conventional_creation=True)

    # This will encode the message of length 10, "a         ".
    input_msg = [ord('a'), 32, 32, 32, 32, 32, 32, 32, 32, 32]
    print(f"Original message: {input_msg}")
    encoded_mesg = grs.encode(input_msg)
    print(f"Encoded message: {encoded_mesg}\n")

    # First consider the case where no error occurred.
    print("No Error:")
    print(f"Syndrome of received message with out error: {grs.syndrome(encoded_mesg)}")
    print(f"Decoded message without error: {grs.decode(encoded_mesg)}\n")

    # By changing out the following lines with another we can see how 1 vs 2 vs 3 errors affects the result
    # Remember for this case we can correct a max of 2 errors
    error = [107, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # One Error
    # error = [0, 4, 0, 250, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Two Errors
    # error = [0, 3, 0, 0, 0, 2, 0, 0, 0, 5, 0, 0, 0, 0, 0]  # Three Errors

    msg_w_error = []
    for m in range(len(encoded_mesg)):
        msg_w_error.append(field.add(encoded_mesg[m], error[m]))

    print("Error:")
    print(f"Received message with error: {msg_w_error},\n\tError: {error}\n")

    # The syndrome of both the msg and the error alone
    print(f"Syndrome of received message with error: {grs.syndrome(msg_w_error)}")
    print(f"Syndrome of error alone: {grs.syndrome(error)}\n")

    decoded_msg = grs.decode(msg_w_error)
    print(f"Decoded received message with error: {decoded_msg}\n")

    print(f"Success: {input_msg == decoded_msg}")





