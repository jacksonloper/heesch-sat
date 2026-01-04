import { useState, useEffect } from 'react'
import WitnessList from './components/WitnessList'
import WitnessViewer from './components/WitnessViewer'
import './App.css'

function App() {
  const [witnesses, setWitnesses] = useState([])
  const [selected, setSelected] = useState(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(null)

  useEffect(() => {
    // Load from static JSONL file built at build time
    fetch('/data/witnesses.jsonl')
      .then(res => {
        if (!res.ok) throw new Error(`Failed to load witnesses: ${res.status}`)
        return res.text()
      })
      .then(text => {
        const polyforms = text
          .trim()
          .split('\n')
          .filter(line => line.length > 0)
          .map(line => JSON.parse(line))
        setWitnesses(polyforms)
        setLoading(false)
      })
      .catch(err => {
        setError(err.message)
        setLoading(false)
      })
  }, [])

  if (loading) {
    return <div className="app loading">Loading witnesses...</div>
  }

  if (error) {
    return <div className="app error">Error: {error}</div>
  }

  return (
    <div className="app">
      <header>
        <div className="header-content">
          <div>
            <h1>Heesch Witness Gallery</h1>
            <p className="subtitle">
              A gallery of interesting polyform tilings and their Heesch numbers
            </p>
          </div>
        </div>
      </header>

      <main>
        <WitnessList
          witnesses={witnesses}
          selected={selected}
          onSelect={setSelected}
        />

        {selected && (
          <WitnessViewer
            witness={selected}
            onClose={() => setSelected(null)}
          />
        )}
      </main>

      <footer className="attribution">
        <p className="gallery-credit">
          This gallery was produced by Jackson Loper to visualize examples from the work of Craig Kaplan.
        </p>
        <h3>References</h3>
        <ul>
          <li>
            <a href="https://github.com/isohedral/heesch-sat" target="_blank" rel="noopener noreferrer">heesch-sat</a> — Source code for computing Heesch numbers of unmarked polyforms using a SAT solver
          </li>
          <li>
            Craig S. Kaplan. Heesch numbers of unmarked polyforms. <em>Contributions to Discrete Mathematics</em>, 17(2):150–171, 2022.{' '}
            <a href="https://cdm.ucalgary.ca/article/view/72886" target="_blank" rel="noopener noreferrer">Available online</a>
          </li>
          <li>
            Craig S. Kaplan. Detecting isohedral polyforms with a SAT solver. <em>GASCom 2024 Abstracts</em>, 118–122, 2024.{' '}
            <a href="https://cgi.cse.unsw.edu.au/~eptcs/paper.cgi?GASCom2024.25" target="_blank" rel="noopener noreferrer">Available online</a>
          </li>
          <li>
            Joseph Myers&apos;s web page about tiling properties of polyforms:{' '}
            <a href="https://www.polyomino.org.uk/mathematics/polyform-tiling/" target="_blank" rel="noopener noreferrer">polyomino.org.uk</a>
          </li>
        </ul>
      </footer>
    </div>
  )
}

export default App
