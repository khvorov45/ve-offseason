# Small (<500) samples

Low (<5%) disease probability and a small (<500) sample size will result in the generation of samples that do not contain any (un)vaccinated cases.

Conditioning on there being at least one vaccinated and at least one unvaccinated case leads to bias because the most likely outcome (1 of each) gives unrepresentative odds in the cases.

Conditioning on there being at least one unvaccinated case leads to bias because all the instances of >0 vaccinated cases but 0 unvaccinated cases are thrown out.

Conclusion: don't do very small samples.
