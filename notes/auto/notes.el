(TeX-add-style-hook
 "notes"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "a4paper" "12pt")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art12"
    "amsmath"
    "times")
   (TeX-add-symbols
    '("ppoly" 3)
    "D")
   (LaTeX-add-labels
    "equ:koornwinder:1"
    "equ:jacobi:1"
    "equ:legendre"
    "equ:sphere:1"))
 :latex)

