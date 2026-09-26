"""Fail-closed NumPy evaluator for a tested subset of MISS_HIT MATLAB ASTs.
Not a MATLAB runtime: unsupported syntax/builtins raise, never silently skip.
"""
from pathlib import Path
import operator, hashlib
import numpy as np
from miss_hit_core.config import Config
from miss_hit_core.errors import Message_Handler
from miss_hit_core.m_language import MATLAB_Latest_Language
from miss_hit_core.m_lexer import MATLAB_Lexer, Token_Buffer
from miss_hit_core.m_parser import MATLAB_Parser


class Return(Exception):
    pass


class Continue(Exception):
    pass


class Break(Exception):
    pass


def scalar(x):
    return np.asarray(x).item()


def integer(x):
    return int(scalar(x))


def mat(x):
    return np.atleast_2d(x)


def truth(x):
    a = np.asarray(x)
    assert a.size == 1, ("non-scalar condition", a.shape)
    return bool(a.item())


def axis(x):
    return next((i for i, n in enumerate(mat(x).shape) if n != 1), 0)


def parse(path):
    path = Path(path)
    mh = Message_Handler("debug")
    mh.register_file(str(path))
    cfg = Config()
    lang = MATLAB_Latest_Language()
    buf = Token_Buffer(MATLAB_Lexer(lang, mh, path.read_text(), str(path)), cfg)
    return MATLAB_Parser(mh, buf, cfg).parse_file()


def matlab_error(identifier, message):
    raise ValueError(f"{identifier}: {message}")


class Engine:
    def __init__(self):
        self.functions = {}
        self.builtins = {}
        self.hashes = {}
        self.calls = {}
        self.builtins.update(
            {
                "size": lambda x, d: mat(x).shape[integer(d) - 1],
                "length": lambda x: 0 if np.size(x) == 0 else max(mat(x).shape),
                "numel": lambda x: np.size(x),
                "isempty": lambda x: np.size(x) == 0,
                "isscalar": lambda x: np.size(x) == 1,
                "zeros": lambda *s: np.zeros(tuple(map(integer, s if len(s) > 1 else s * 2))),
                "ones": lambda *s: np.ones(tuple(map(integer, s if len(s) > 1 else s * 2))),
                "eye": lambda n: np.eye(integer(n)),
                "reshape": lambda x, *s: np.asarray(x).reshape(tuple(map(integer, s)), order="F"),
                "repmat": lambda x, *s: np.tile(mat(x), tuple(map(integer, s))),
                "sum": lambda x, d=None: np.sum(
                    mat(x), axis=axis(x) if d is None else integer(d) - 1, keepdims=True
                ),
                "any": lambda x: np.any(mat(x), axis=axis(x), keepdims=True),
                "all": lambda x: np.all(mat(x), axis=axis(x), keepdims=True),
                "cumsum": lambda x: np.cumsum(mat(x), axis=axis(x)),
                "sqrt": np.sqrt,
                "cos": np.cos,
                "acos": np.arccos,
                "isnan": np.isnan,
                "isfinite": np.isfinite,
                "error": matlab_error,
                "norm": lambda x: np.linalg.norm(x),
                "trace": np.trace,
                "floor": np.floor,
                "ceil": np.ceil,
                "real": np.real,
                "isreal": np.isrealobj,
                "gradient": lambda x: np.gradient(x, axis=1),
                "find": lambda x: (np.flatnonzero(np.asarray(x).ravel(order="F")) + 1).reshape(
                    -1, 1
                ),
                "strcmp": lambda a, b: a == b,
                "unique": lambda x: np.unique(x).reshape(-1, 1),
                "randperm": lambda n: (np.arange(integer(n)) + 1).reshape(1, -1),
                "fprintf": lambda *a: None,
                "disp": lambda *a: None,
                "isfield": lambda s, k: k in s,
                "sprintf": lambda fmt, *args: fmt
                % tuple(scalar(x) if isinstance(x, np.ndarray) else x for x in args),
                "fullfile": lambda *parts: str(Path(*parts)),
            }
        )

    def load(self, path):
        p = Path(path)
        tree = parse(p)
        self.hashes[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
        for f in tree.l_functions:
            self.functions[str(f.n_sig.n_name)] = f

    def call(self, name, args, nout=1):
        self.calls[name] = self.calls.get(name, 0) + 1
        if name in self.builtins:
            return self.builtins[name](*args)
        f = self.functions[name]
        env = {
            "nargin": len(args),
            "nargout": nout,
            "eps": np.finfo(float).eps,
            "NaN": np.nan,
            "Inf": np.inf,
        }
        env.update(
            {
                str(n): x.copy() if isinstance(x, np.ndarray) else x
                for n, x in zip(f.n_sig.l_inputs, args)
            }
        )
        try:
            self.block(f.n_body, env)
        except Return:
            pass
        outs = [env[str(x)] for x in f.n_sig.l_outputs[:nout]]
        return outs[0] if nout == 1 else tuple(outs)

    def indices(self, arr, nodes, env):
        shape = mat(arr).shape
        values = []
        for i, n in enumerate(nodes):
            size = np.size(arr) if len(nodes) == 1 else shape[i]
            if type(n).__name__ == "Reshape":
                v = np.arange(size)
            else:
                e = dict(env, end=size)
                v = np.asarray(self.expr(n, e))
                v = (
                    np.flatnonzero(v.ravel(order="F"))
                    if v.dtype == bool
                    else v.astype(int).ravel(order="F") - 1
                )
            if np.any(v < 0) or np.any(v >= size):
                raise IndexError((v, size))
            values.append(v)
        return values

    def read(self, arr, nodes, env):
        idx = self.indices(arr, nodes, env)
        a = np.asarray(arr)
        if len(idx) == 1:
            value = a.ravel(order="F")[idx[0]]
            return (
                value.reshape(1, -1)
                if a.ndim == 2 and a.shape[0] == 1 and type(nodes[0]).__name__ != "Reshape"
                else value.reshape(-1, 1)
            )
        return mat(a)[np.ix_(*idx)]

    def assign(self, n, val, e):
        kind = type(n).__name__
        if kind == "Identifier":
            if str(n) != "~":
                e[str(n)] = val.copy() if isinstance(val, np.ndarray) else val
        elif kind == "Reference":
            name = str(n.n_ident)
            a = mat(e[name]).copy()
            idx = self.indices(a, n.l_args, e)
            if np.size(val) == 0:
                assert len(idx) == 2 and type(n.l_args[0]).__name__ == "Reshape"
                e[name] = np.delete(a, idx[1], axis=1)
                return
            if len(idx) == 1:
                flat = a.ravel(order="F").copy()
                flat[idx[0]] = np.asarray(val).ravel(order="F")
                e[name] = flat.reshape(a.shape, order="F")
            else:
                a[np.ix_(*idx)] = val
                e[name] = a
        else:
            raise NotImplementedError(("assignment", kind))

    def expr(self, n, e, nout=1):
        k = type(n).__name__
        if k == "Identifier":
            return e[str(n)]
        if k == "Number_Literal":
            return float(str(n))
        if k in ("Char_Array_Literal", "String_Literal"):
            return n.t_string.value
        if k in ("Selection", "Dynamic_Selection"):
            return self.expr(n.n_prefix, e)[
                str(n.n_field) if k == "Selection" else self.expr(n.n_field, e)
            ]
        if k in ("Matrix_Expression", "Cell_Expression"):
            rows = [[self.expr(x, e) for x in row.l_items] for row in n.n_content.l_items]
            if k == "Cell_Expression":
                return [x for row in rows for x in row]
            if not rows or not any(rows):
                return np.empty((0, 0))
            if all(isinstance(x, str) for row in rows for x in row):
                return "".join(x for row in rows for x in row)
            return np.vstack([np.hstack([mat(x) for x in row]) for row in rows if row])
        if k == "Range_Expression":
            first = scalar(self.expr(n.n_first, e))
            last = scalar(self.expr(n.n_last, e))
            step = scalar(self.expr(n.n_stride, e)) if n.n_stride else 1
            return np.arange(first, last + step * 1e-10, step).reshape(1, -1)
        if k == "Cell_Reference":
            value = self.read(self.expr(n.n_ident, e), n.l_args, e)
            if value.size != 1:
                raise NotImplementedError("cell comma-separated list")
            return value.item()
        if k == "Reference" and type(n.n_ident).__name__ != "Identifier":
            return self.read(self.expr(n.n_ident, e), n.l_args, e)
        if k in ("Reference", "Function_Call"):
            name = str(n.n_ident if k == "Reference" else n.n_name)
            if name in e:
                return self.read(e[name], n.l_args, e)
            return self.call(name, [self.expr(x, e) for x in n.l_args], nout)
        if k == "Unary_Operation":
            x = self.expr(n.n_expr, e)
            op = n.t_op.value
            return {
                "'": lambda: mat(x).conj().T,
                ".'": lambda: mat(x).T,
                "-": lambda: -x,
                "+": lambda: x,
                "~": lambda: np.logical_not(x),
            }[op]()
        if k in ("Binary_Operation", "Binary_Logical_Operation"):
            op = n.t_op.value
            a = self.expr(n.n_lhs, e)
            if op == "||":
                return truth(a) or truth(self.expr(n.n_rhs, e))
            if op == "&&":
                return truth(a) and truth(self.expr(n.n_rhs, e))
            b = self.expr(n.n_rhs, e)
            if op == "*":
                return a * b if np.size(a) == 1 or np.size(b) == 1 else mat(a) @ mat(b)
            if op == "/":
                if np.size(b) == 1:
                    return a / b
                a, b = mat(a), mat(b)
                if b.shape[0] != b.shape[1]:
                    raise NotImplementedError("only square nonsingular right division is supported")
                return np.linalg.solve(b.T, a.T).T
            ops = {
                "+": operator.add,
                "-": operator.sub,
                ".*": operator.mul,
                "./": operator.truediv,
                ".^": operator.pow,
                "^": operator.pow,
                "==": operator.eq,
                "~=": operator.ne,
                "<": operator.lt,
                ">": operator.gt,
                "<=": operator.le,
                ">=": operator.ge,
                "&": np.logical_and,
                "|": np.logical_or,
            }
            return ops[op](a, b)
        raise NotImplementedError(("expression", k, str(n)))

    def block(self, b, e):
        for n in b.l_statements:
            k = type(n).__name__
            try:
                if k == "Simple_Assignment_Statement":
                    self.assign(n.n_lhs, self.expr(n.n_rhs, e), e)
                elif k == "Compound_Assignment_Statement":
                    values = self.expr(n.n_rhs, e, len(n.l_lhs))
                    for dst, val in zip(n.l_lhs, values):
                        self.assign(dst, val, e)
                elif k == "General_For_Statement":
                    for x in np.asarray(self.expr(n.n_expr, e)).ravel(order="F"):
                        e[str(n.n_ident)] = x
                        try:
                            self.block(n.n_body, e)
                        except Continue:
                            continue
                        except Break:
                            break
                elif k == "If_Statement":
                    for a in n.l_actions:
                        if a.n_expr is None or truth(self.expr(a.n_expr, e)):
                            self.block(a.n_body, e)
                            break
                elif k == "Switch_Statement":
                    val = self.expr(n.n_expr, e)
                    for a in n.l_actions:
                        cases = self.expr(a.n_expr, e) if a.n_expr else None
                        if cases is None or val in (cases if isinstance(cases, list) else [cases]):
                            self.block(a.n_body, e)
                            break
                elif k == "Naked_Expression_Statement":
                    x = n.n_expr
                    if type(x).__name__ == "Function_Call" and str(x.n_name) in (
                        "clear",
                        "save",
                        "load",
                    ):
                        args = [self.expr(v, e) for v in x.l_args]
                        if str(x.n_name) == "clear":
                            for key in args:
                                e.pop(key, None)
                        elif str(x.n_name) == "load":
                            e.update(self.builtins["load_workspace"](args))
                        else:
                            self.builtins["save_workspace"](args, e)
                    else:
                        self.expr(x, e)
                elif k == "Return_Statement":
                    raise Return()
                elif k == "Continue_Statement":
                    raise Continue()
                elif k == "Break_Statement":
                    raise Break()
                else:
                    raise NotImplementedError(("statement", k))
            except (Return, Continue, Break):
                raise
            except Exception as ex:
                ex.add_note(str(n.loc()))
                raise
