# Contributing to MamPipe

Thank you for your interest in contributing to MamPipe!

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue with:
- A clear description of the problem
- Steps to reproduce
- Expected vs actual behavior
- Your environment (OS, Nextflow version, etc.)

### Suggesting Enhancements

We welcome suggestions for new features or improvements. Please open an issue describing:
- The enhancement you'd like to see
- Why it would be useful
- Any implementation ideas you have

### Pull Requests

1. Fork the repository
2. Create a new branch (`git checkout -b feature/your-feature`)
3. Make your changes
4. Test your changes thoroughly
5. Commit with clear messages (`git commit -m 'Add feature X'`)
6. Push to your fork (`git push origin feature/your-feature`)
7. Open a Pull Request

### Code Style

- Follow existing code formatting
- Comment complex logic
- Update documentation for new features
- Add appropriate error handling

### Testing

Before submitting a PR:
- Test locally with sample data
- Ensure all processes complete successfully
- Check that output files are generated correctly

## Development Setup

```bash
# Clone your fork
git clone https://github.com/YOUR_USERNAME/mampipe.git
cd mampipe

# Create a test branch
git checkout -b test-feature

# Make changes and test
nextflow run main.nf --reads test_data/*.fastq --ref test_data/*.fna
```

## Questions?

Feel free to open an issue for any questions about contributing.
