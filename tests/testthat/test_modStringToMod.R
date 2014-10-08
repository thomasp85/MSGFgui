context('Parsing of modification strings from client')

modstring1 <- 'N:Carbamidomethyl;C:C2H3N1O1;R:C;T:fix;P:any'
modstring2 <- 'N:Oxidation;W:15.994915;R:M;T:opt;P:any'
mod1 <- MSGFplus::msgfParModification(name='Carbamidomethyl', composition='C2H3N1O1', residues='C', type='fix', position='any')
mod2 <- MSGFplus::msgfParModification(name='Oxidation', mass=15.994915, residues='M', type='opt', position='any')

test_that('parsing works', {
    expect_equal(modStringToMod(modstring1), mod1)
    expect_equal(modStringToMod(modstring2), mod2)
})